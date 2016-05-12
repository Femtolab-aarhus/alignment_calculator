#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#   Copyright 2016 Anders Aspegren Søndergaard / Femtolab, Aarhus University
#
#   This file is part of Alignment calculator.
#
#   Alignment calculator is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Alignment calculator is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with Alignment calculator. If not, see <http://www.gnu.org/licenses/>.


import numpy
import numpy as np
import scipy
import scipy.integrate
import scipy.sparse
from scipy.integrate import ode
from scipy import linalg
from numpy import exp,log,sqrt,pi
import interaction
import U2dcalc
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
can_propagate_using_ODE = False;
try: # to load the C library for propagation.
    if (utils.running_on_windows()):
        libpropagation = ctypes.CDLL("./libpropagation.dll");
    else:
        libpropagation = ctypes.CDLL("./libpropagation.so");
    cplx_ndptr = ndpointer(flags=("CONTIGUOUS",),dtype=numpy.complex);
    real_ndptr = ndpointer(flags=("CONTIGUOUS",),dtype=numpy.double);
    libpropagation.fieldfree_propagation.restype = None;
    libpropagation.fieldfree_propagation.argtypes = (ct.c_int, cplx_ndptr, ct.c_double, ct.c_size_t, real_ndptr, real_ndptr, real_ndptr, real_ndptr, real_ndptr, cplx_ndptr, real_ndptr, ct.c_bool, real_ndptr, real_ndptr);
    libpropagation.propagate_field.restype = ct.c_int;
    libpropagation.propagate_field.argtypes = (ct.c_size_t, ct.c_size_t, ct.c_size_t, ct.c_double, ct.c_double, ct.c_double, ct.c_double, cplx_ndptr, real_ndptr, real_ndptr, cplx_ndptr, cplx_ndptr);
    
    # Don't use the ODE method if it is not available, i.e. if we compiled
    # without GSL dependence.
    try:
        libpropagation.propagate_field_ODE.restype = ct.c_int;
        libpropagation.propagate_field_ODE.argtypes = (ct.c_size_t, ct.c_size_t, ct.c_double, ct.c_double, ct.c_double, ct.c_double, cplx_ndptr, real_ndptr, real_ndptr, real_ndptr, real_ndptr, ct.c_double, ct.c_double);
        can_propagate_using_ODE = True;
    except:
        pass;

except OSError: # Fall back to the python implementation
    print("Not using the C propagation library (compile it with make). Propagation will be slow.",file=sys.stderr);


def transfer_JKM(J,K,M,KMsign,Jmax,peak_intensity,FWHM,t,molecule,use_ODE=True):
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
     use_ODE: Use the ODE solver instead of the matrix method. This is typically
              fastest, particularly for short pulses.
     '''

     psi = numpy.zeros(Jmax+1,dtype=numpy.complex);
     psi[J] = 1;

     return transfer_KM(psi,K,M,KMsign,Jmax,peak_intensity,FWHM,t,molecule,use_ODE);



def transfer_KM(psi,K,M,KMsign,Jmax,peak_intensity,FWHM,t,molecule,use_ODE=True,abstol=1e-8,reltol=1e-8):
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
     use_ODE: Use the ODE solver instead of the matrix method. This is typically
              fastest, particularly for short pulses.
     abstol,reltol: absolute and relative total error tolerances for the
                    propagation over the pulse. The larger the tolerance
                    for error, the faster the execution time. 
                    Only used when use_ODE is True.
     '''
     A = molecule.A;
     B = molecule.B;
     alpha_perp = molecule.alpha_perp_volume;
     delta_alpha = molecule.alpha_par_volume-alpha_perp;
     
     hbar=1.05457173e-34;
     c = 299792458.0;

     sigma = FWHM/(2*sqrt(2*log(2)));

     Js = numpy.arange(0.0,Jmax+1,1);
     E_rot = B*Js*(Js+1) + (A-B)*K**2; # *hbar / hbar for solver

     # Note: since V and E_0_max ends up getting multiplied together,
     # 2 factors of epsilon_0 cancels
     # (since alpha is given in polarizability volume in Angstrom^3 =
     # 1/4pi epsilon_0 times polarizability )
     # Therefore, these factors are skipped.
    
     E_0_squared_max = 2*peak_intensity/c;

     if (not can_propagate_using_ODE):
         use_ODE = False;

     if (not use_ODE or Jmax<=2):
        V0, V1, V2, eig, vec = interaction.interaction_matrix_div_E0_squared(Jmax,K,M,KMsign,delta_alpha,alpha_perp,diagonalize=True);
        # Divide by hbar for solver
        eig = eig*4*pi*1e-30/hbar
        psi_t = propagate(psi,t,E_rot,E_0_squared_max,sigma,eig,vec);
     else:
        V0, V1, V2, = interaction.interaction_matrix_div_E0_squared(Jmax,K,M,KMsign,delta_alpha,alpha_perp,diagonalize=False);
        V0 = V0*4*pi*1e-30/hbar; # Divide by hbar for solver
        V1 = V1*4*pi*1e-30/hbar;
        V2 = V2*4*pi*1e-30/hbar;
        psi_t = propagate_ODE(psi,t,E_rot,E_0_squared_max,sigma,V0,V1,V2,abstol,reltol);

     return psi_t;


def propagate(psi_0,time,E_rot,E_0_squared_max,sigma,eig,vec):
     ''' Propagate the wave function psi_0 (given as coefficients)
     through a gaussian laser pulse centered at time t = 0.
     
     psi_0: initial wave function
     time: array of timesteps to integrate. psi_0 is given at time[0].
           The time step must be constant!
     E_rot: array of rotational energies for each J.
     E_0_squared_max: Max value of the square of the E field during pulse
     sigma: temporal variance of the gaussian
     eig, vec: diagonalized interaction matrix V = vec * diag(eig) * vec.T
     
     '''
     if (numpy.any(numpy.abs(numpy.diff(numpy.diff(time)))>1e-20)):
         raise RuntimeError("Pulse time steps must be equidistant.");

     # Within the maximum time step, a Gaussian is approximately constant.
     dt_max = 2.4*sigma/150; 
     dt = time[1]-time[0];
     # If our given time step is larger than this, use a smaller time step
     if (dt > dt_max):
         scale = int(numpy.ceil(dt/dt_max));
         dt = dt/scale;
     else:
         scale = 1;
    
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
                if (k > 0): # Double half step:
                    psi_t[i,:] = expRot*psi_t[i,:];
                tp = t+k*dt;
                E_0_squared = E_0_squared_max * numpy.exp(-(tp**2)/(2*sigma**2)); 
                psi_t[i,:] = ((numpy.exp(-dt*E_0_squared*1j*eig))*vec).dot(vec.T.dot(psi_t[i,:]));
                if (numpy.max(numpy.abs(psi_t[i,-2:])**2) > 1e-5):
                    raise RuntimeError("Basis size too small");
            psi_t[i,:] = expRot2*psi_t[i,:];
            i = i + 1;
 
     return psi_t;

def propagate_ODE(psi_0,time,E_rot,E_0_squared_max,sigma,V0,V1,V2,abstol=1e-8,reltol=1e-8):
     ''' Same as propagate, except it uses an ODE solver instead
     of the matrix method.
     
     psi_0: initial wave function
     time: array of timesteps to integrate. psi_0 is given at time[0].
           The time step must be constant!
     E_rot: array of rotational energies for each J.
     E_0_squared_max: Max value of the square of the E field during pulse
     sigma: temporal variance of the gaussian
     V0,V1,V2: the three bands of the symmetric 5 diagonal interaction matrix,
               V0 being the diagonal.
     abstol, reltol: Error tolerances used for the ODE solver.
     
     '''
     if (numpy.any(numpy.abs(numpy.diff(numpy.diff(time)))>1e-20)):
         raise RuntimeError("Pulse time steps must be equidistant.");
     dt = time[1]-time[0];

     psi_t = numpy.empty((len(time),len(psi_0)), dtype=numpy.complex)
     psi_t[0,:] = psi_0;
     try:
         res = libpropagation.propagate_field_ODE(len(time),len(psi_0),time[0],dt,E_0_squared_max,sigma,psi_t,V0,V1,V2,E_rot, abstol, reltol);
     except:
         raise RuntimeError("For ODE propagation, you need the libpropagation C library, compiled with GSL support.");
     if (res != 0):
         raise RuntimeError("Basis size too small");
 
     return psi_t;


def fieldfree_propagation(psi_0,t0,times,E_rot,Jmax,K,M,KMsign,do_cos2d=False):

     U, Udiag, U1, U2 = interaction.MeanCos2Matrix(Jmax,K,M,KMsign);
     cos2d = numpy.array([]);
     if (do_cos2d):
        U2d = U2dcalc.MeanCos2dMatrix(Jmax,K,M,KMsign);
     else:
        U2d = U;
    
     if (libpropagation and Jmax >= 3): # Call a C function instead of using numpy.
        psi = numpy.empty((len(times),Jmax+1),dtype=numpy.complex);
        cos2 = numpy.empty((len(times),),dtype=numpy.double);
        if (do_cos2d):
            cos2d = numpy.empty((len(times),),dtype=numpy.double);
        libpropagation.fieldfree_propagation(Jmax,psi_0,t0,len(times),times,E_rot,Udiag,U1,U2,psi,cos2,do_cos2d,U2d,cos2d);
        return psi,cos2,cos2d;
    
     U = scipy.sparse.diags([Udiag,U1,U1,U2,U2],[0,-1,1,-2,2],[len(Udiag)]*2)

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

