#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# This file is kept as reference. It contains the propagate_old method
# which propagates the wave function by solving an ODE for the 
# expansion coefficients. It was slow and tended to be inaccurate for long
# pulses, but it is kept because the method itself is interesting.
# If implemented in C, it could perhaps even beat the split step method.
# It too, suffered from repeated memory allocation and deallocation
# in the jacobian method which gets called in a tight loop.

import numpy
import numpy as np
import scipy
import scipy.integrate
import scipy.sparse
from scipy.integrate import ode
from scipy import linalg
from numpy import exp,log,sqrt,pi
import interaction
import datetime
import time
import io
import multiprocessing
import multiprocessing.sharedctypes
import ctypes
import ctypes as ct
from numpy.ctypeslib import ndpointer
import sys
from itertools import repeat
import tempfile
import utils

libpropagation = None;
try: # to load the C library for fieldfree propagation.
    if (utils.running_on_windows()):
        libpropagation = ctypes.CDLL("./libpropagation.dll");
    else:
        libpropagation = ctypes.CDLL("./libpropagation.so");
    cplx_ndptr = ndpointer(flags="CONTIGUOUS",dtype=numpy.complex);
    real_ndptr = ndpointer(flags="CONTIGUOUS",dtype=numpy.double);
    libpropagation.fieldfree_propagation.restype = None;
    libpropagation.fieldfree_propagation.argtypes = (ct.c_int, cplx_ndptr, ct.c_double, ct.c_size_t, real_ndptr, real_ndptr, real_ndptr, real_ndptr, real_ndptr, cplx_ndptr, real_ndptr, ct.c_bool, real_ndptr, real_ndptr);
    libpropagation.propagate_field.restype = ct.c_int;
    libpropagation.propagate_field.argtypes = (ct.c_size_t, ct.c_size_t, ct.c_size_t, ct.c_double, ct.c_double, ct.c_double, ct.c_double, cplx_ndptr, real_ndptr, real_ndptr, cplx_ndptr, cplx_ndptr);
except OSError: # Fall back to the python implementation
    print("Not using the C library (compile it with make).",file=sys.stderr);


def transfer_JKM(J,K,M,KMsign,Jmax,peak_intensity,FWHM,t,molecule):
     ''' Calculate the wave function after a laser pulse acting on a molecule
     in the state |JKM>.

     J, K, M: J, K and M quantum numbers (must be nonnegative)
     KMsign: Relative sign between K and M.
     Jmax: Maximum J state to include
     peak_intensity: Peak intensity of laser pulse
     FWHM: Full width at half maximum of laser pulse
     t: Solve Schrödinger equation from from t[0] to t[1],t[2],...,t[-1]
             (the pulse is centered at t=0)
     molecule: Configuration object for the molecule.
     '''

     psi = numpy.zeros(Jmax+1,dtype=numpy.complex);
     psi[J] = 1;

     return transfer_KM(psi,K,M,KMsign,Jmax,peak_intensity,FWHM,t,molecule);


def transfer_KM(psi,K,M,KMsign,Jmax,peak_intensity,FWHM,t,molecule):
     ''' Calculate the wave function after a laser pulse acting on a molecule
     in the state \sum_J C_J|KM>.

     psi: initial state wavefunction (C_J coefficients)
     K, M: K and M quantum numbers (must be nonnegative)
     KMsign: Relative sign between K and M.
     Jmax: Maximum J state to include
     peak_intensity: Peak intensity of laser pulse
     FWHM: Full width at half maximum of laser pulse
     t: Solve Schrödinger equation from from t[0] to t[1],t[2],...,t[-1]
             (the pulse is centered at t=0)
     molecule: Configuration object for the molecule.
     '''
     A = molecule.A;
     B = molecule.B;
     alpha_perp = molecule.alpha_perp_volume;
     delta_alpha = molecule.alpha_par_volume-alpha_perp;
     
     hbar=1.05457173e-34;
     c = 299792458.0;

     sigma = FWHM/(2*sqrt(2*log(2)));

     # Note: since V and E_0_max ends up getting multiplied together,
     # 2 factors of epsilon_0 cancels
     # (since alpha is given in polarizability volume in Angstrom^3 =
     # 1/4pi epsilon_0 times polarizability )
     # Therefore, these factors are skipped.
     V0, V1, V2 = [mat*4*pi*1e-30 for mat in interaction.interaction_matrix_div_E0_squared(Jmax,K,M,KMsign,delta_alpha,alpha_perp)];
     E_0_squared_max = 2*peak_intensity/c;
     
     V0 /= hbar; # Divide by hbar for solver
     V1 /= hbar;
     V2 /= hbar;

     Js = numpy.arange(0,Jmax+1,1);
     E_rot = B*Js*(Js+1) + (A-B)*K**2; # *hbar / hbar for solver
     
     psi_t = propagate(psi,t,E_rot,E_0_squared_max,sigma,V0,V1,V2);
     #psi_t = propagate_old(psi,t,E_rot,E_0_squared_max,sigma,V0,V1,V2);

     return psi_t;


def propagate(psi_0,time,E_rot,E_0_squared_max,sigma,V0,V1,V2,tol=None):
     ''' Propagate the wave function psi_0 (given as coefficients)
     through a gaussian laser pulse centered at time t = 0.
     
     psi_0: initial wave function
     time: array of timesteps to integrate. psi_0 is given at time[0].
           The time step must be constant!
     E_rot: array of rotational energies for each J.
     E_0_squared_max: Max value of the square of the E field during pulse
     sigma: temporal variance of the gaussian
     V: Interaction matrix (symmetric 5 diagonal)
     
     '''
     if (numpy.any(numpy.abs(numpy.diff(numpy.diff(time)))>1e-20)):
         raise RuntimeError("Pulse time steps must be equidistant.");

     # Within the minimum time step, a Gaussian is approximately constant.
     dt_min = 2.4*sigma/150; 
     dt = time[1]-time[0];
     # If our given time step is larger than this, use a smaller time step
     if (dt > dt_min):
         scale = int(numpy.ceil(dt/dt_min));
         dt = dt/scale;
     else:
         scale = 1;
     cat = numpy.concatenate;

     a = numpy.vstack((cat(([0,0],V2)),cat(([0],V1)),V0));
    
     eig,vec = linalg.eig_banded(a);
     vec = numpy.ascontiguousarray(vec);
    
     expRot = numpy.exp(-1j*dt*E_rot);
     expRot2 = numpy.exp(-1j*dt*E_rot/2);
    
     psi_t = numpy.empty((len(time),len(psi_0)), dtype=numpy.complex)
     psi_t[0,:] = psi_0;
     if (libpropagation):
         res = libpropagation.propagate_field(len(time),scale,len(psi_0),time[0],dt,E_0_squared_max,sigma,psi_t,eig,vec,expRot,expRot2);
         if (res != 0):
             raise RuntimeError("Basis size too small");
     else:
         i = 1;
         for t in time[:-1]: # use split step
            psi_t[i,:] = expRot2*psi_t[i-1,:];
            for k in range(scale):
                # I(t), gaussian beam
                if (k > 0): # Double step:
                    psi_t[i,:] = expRot*psi_t[i,:];
                tp = t+k*dt;
                E_0_squared = E_0_squared_max * numpy.exp(-(tp**2)/(2*sigma**2)); 
                psi_t[i,:] = ((numpy.exp(-dt*E_0_squared*1j*eig))*vec).dot(vec.T.dot(psi_t[i,:]));
                if (numpy.max(numpy.abs(psi_t[i,-2:])**2) > 1e-5):
                    raise RuntimeError("Basis size too small");
            psi_t[i,:] = expRot2*psi_t[i,:];
            i = i + 1;
 
     return psi_t;


def propagate_old(psi_0,time,E_rot,E_0_squared_max,sigma,V0,V1,V2,tol=None):
     ''' Propagate the wave function psi_0 (given as coefficients)
     through a gaussian laser pulse centered at time t = 0.
     
     psi_0: initial wave function
     time: array of timesteps to integrate. psi_0 is given at time[0].
     E_rot: array of rotational energies for each J.
     E_0_squared_max: Max value of the square of the E field during pulse
     sigma: temporal variance of the gaussian
     V: Interaction matrix (symmetric 5 diagonal)
     
     '''
     def ddt(t,y):
          return dPsi_dt(t,y,E_0_squared_max, E_rot, sigma, V0, V1, V2);
     def D(t,y):
          return jac(t,y,E_0_squared_max, E_rot, sigma, V0, V1, V2);
     
     if (tol is None):
          Jmax = V0.shape[0]-1;
          tol = 1e-5/(Jmax**2);
     ode = scipy.integrate.complex_ode(ddt,jac=D);
     #ode = scipy.integrate.ode(dPsi_dt,jac=jac);
     #ode.set_integrator('vode', method='bdf',first_step=dt,nsteps=2**31-1,rtol=tol,atol=tol);
     ode.set_integrator('vode', method='bdf',nsteps=2**31-1,rtol=tol,atol=tol)#,atol=tol);
     #ode.set_integrator('zvode', method='bdf',first_step=dt,nsteps=2**31-1,rtol=tol,atol=tol);
     #ode.set_f_params(E_0_squared_max, E_rot, sigma, V0, V1, V2);
     #ode.set_jac_params(E_0_squared_max, E_rot, sigma, V);
     ode.set_initial_value(psi_0,time[0]);
     
     psi_t = numpy.empty((len(time),len(psi_0)), dtype=numpy.complex)
     psi_t[0,:] = psi_0;
     i = 1;
     for t in time[1:]:
          psi_t[i,:] = ode.integrate(t);
          i+=1;
          if (not ode.successful()):
               raise RuntimeError("Failed to integrate to t = "+str(t)+".");

          #print ode.t#, abs(psi_t[0:5])**2;
     
     if (numpy.max(numpy.abs(psi_t[:,-2:])**2) > 1e-5):
         raise RuntimeError("Basis size too small");

     return psi_t;


def dPsi_dt(t,psi,peak_field_amplitude_squared,E_rot,sigma,V0,V1,V2):
     ''' Interaction with a laser pulse centered at t=0
     
     E_rot: eigenvalues of the rotational Hamiltonian divided by hbar
     
     V: Interaction matrix divided by E_0 (a constant!) and by hbar

     Diff equation autonomous: choose max pulse to occur at t=0.

     '''
     
     E_0_squared = peak_field_amplitude_squared * \
                   exp(-t**2/(2*sigma**2)); # I(t), gaussian beam
     
     #dPsi = numpy.zeros_like(psi);
     
     #Jmax = len(psi)-1;

     # Rotational term (diagonal):
     dPsi = E_rot*psi;
     #V=numpy.diag(V0) + numpy.diag(V1,1)+numpy.diag(V1,-1)+numpy.diag(V2,2)+numpy.diag(V2,-2);
     # Interaction term:
     dPsi += E_0_squared*V0*psi;
     dPsi[:-1] += E_0_squared*V1*psi[1:]; 
     dPsi[1:] += E_0_squared*V1*psi[:-1];
     dPsi[:-2] += E_0_squared*V2*psi[2:];
     dPsi[2:] += E_0_squared*V2*psi[:-2];
     #for J in range(Jmax+1): # The above is a vectorized version of this:
     #     JPmin = max(0,J-2); # J prime min,max
     #     JPmax = min(Jmax,J+2);
     #     for JP in range(JPmin,JPmax+1): # Interaction term
     #          dPsi[J] += E_0_squared * V[J][JP]*psi[JP];

     return -1j*dPsi; 

def jac(t,psi,peak_field_amplitude_squared,E_rot,sigma,V0,V1,V2):
     E_0_squared = peak_field_amplitude_squared * \
                   exp(-t**2/(2*sigma**2)); # I(t), gaussian beam

     D = numpy.diag(-1j*E_rot*V0);
     D += numpy.diag(-1j*E_0_squared*V1,1);
     D += numpy.diag(-1j*E_0_squared*V1,-1);
     D += numpy.diag(-1j*E_0_squared*V2,2);
     D += numpy.diag(-1j*E_0_squared*V2,-2);

     return D;


def fieldfree_propagation(psi_0,t0,times,E_rot,Jmax,K,M,KMsign,do_cos2d=False):

     U, Udiag, U1, U2 = interaction.MeanCos2Matrix(Jmax,K,M,KMsign);
     cos2d = numpy.array([]);
     if (do_cos2d):
        U2d = interaction.MeanCos2dMatrix(Jmax,K,M,KMsign);
     else:
        U2d = U;
    
     if (libpropagation and Jmax >= 3): # Call a C function instead of using numpy.
        psi = numpy.empty((len(times),Jmax+1),dtype=numpy.complex);
        cos2 = numpy.empty((len(times),),dtype=numpy.double);
        if (do_cos2d):
            cos2d = numpy.empty((len(times),),dtype=numpy.double);
        libpropagation.fieldfree_propagation(Jmax,psi_0,t0,len(times),times,E_rot,Udiag,U1,U2,psi,cos2,do_cos2d,U2d,cos2d);
        return psi,cos2,cos2d;

     U = scipy.sparse.dia_matrix(U); # This appears to be somewhat faster

     phase = exp(-1j*numpy.outer(E_rot,(times-t0)));
     psi = psi_0.reshape(Jmax+1,1)*phase;
    
     # We only need the diagonal in psi^H x U x psi.
     # since diag(AxB) = [sum A1j*Bj1, sum A2j*Bj2 ...] =
     # [sum A^Tj1*Bj1, sum A^Tj2*Bj2 ...] = A^T*B summed along the colums
     cos2 = numpy.real(numpy.conj(psi)*U.dot(psi)).sum(0);
     if (do_cos2d):
         cos2d = numpy.real(numpy.conj(psi)*U2d.dot(psi)).sum(0);

     psi=numpy.transpose(psi);
     
     return psi,cos2,cos2d; # psi[t][n] contains n'th component at time index t

