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
import scipy.interpolate
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
import os
import matplotlib.pyplot as plt

libpropagation = None;
can_propagate_using_ODE = False;
try: # to load the C library for propagation.
    libpropagation = ctypes.CDLL("./libpropagation."+utils.library_file_ending());
    cplx_ndptr = ndpointer(flags=("CONTIGUOUS",),dtype=numpy.complex);
    real_ndptr = ndpointer(flags=("CONTIGUOUS",),dtype=numpy.double);
    libpropagation.fieldfree_propagation.restype = None;
    libpropagation.fieldfree_propagation.argtypes = (ct.c_int, ct.c_int, ct.c_size_t, cplx_ndptr, ct.c_double, ct.c_size_t, real_ndptr, real_ndptr, real_ndptr, real_ndptr, real_ndptr, cplx_ndptr, real_ndptr, ct.c_bool, real_ndptr, real_ndptr);
    libpropagation.propagate_field.restype = ct.c_int;
    libpropagation.propagate_field.argtypes = (ct.c_size_t, ct.c_size_t, ct.c_size_t, ct.c_double, ct.c_double, ct.c_double, ct.c_double, cplx_ndptr, real_ndptr, real_ndptr, real_ndptr, cplx_ndptr, cplx_ndptr, real_ndptr, ct.c_bool);

    # Don't use the ODE method if it is not available, i.e. if we compiled
    # without GSL dependence.
    try:
        libpropagation.propagate_field_ODE.restype = ct.c_int;
        libpropagation.propagate_field_ODE.argtypes = (ct.c_size_t, ct.c_size_t, ct.c_double, ct.c_double, ct.c_double, ct.c_double, cplx_ndptr, real_ndptr, real_ndptr, real_ndptr, real_ndptr, ct.c_double, ct.c_double);
        can_propagate_using_ODE = True; #HARDCODED
    except:
        print("ODE method not compiled into the C library.",file=sys.stderr);
    try:
        if (os.environ["ALIGNMENT_CALC_MAY_USE_ODE_SOLVER"] == "NO"):
            can_propagate_using_ODE = False;
    except:
        pass;

except OSError as e: # Fall back to the python implementation
    print(e,file=sys.stderr);
    print("Not using the C propagation library (compile it with make). Propagation will be slow.",file=sys.stderr);


def linterp(t_pulse,y_pulse,t_prob):
    pulse = scipy.interpolate.interp1d(t_pulse, y_pulse);
    p = pulse(t_prob)
    p[t_prob<t_pulse[0]] = 0;
    p[t_prob>t_pulse[-1]] = 0;
    return p;

def plotlinterp(t_new,alpha,t0,y0):
    y_new=numpy.zeros(t_new.size)
    for i in range(len(t_new)):
        y_new[i]=y0[i]+alpha[i]*(t_new[i]-t0[i])
    plt.plot(t_new,y_new)
    plt.show()


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



def transfer_KM(psi,K,M,KMsign,Jmax,peak_intensity,FWHM,t,molecule,xc_filename,use_ODE=True,abstol=1e-8,reltol=1e-8):
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
     D = molecule.D;
     delta_alpha = molecule.delta_alpha;
     alpha_perp = molecule.alpha_perp_volume

     hbar=1.05457173e-34;
     c = 299792458.0;

     sigma = FWHM/(2*sqrt(2*log(2)));

     Js = numpy.arange(0.0,Jmax+1,1);
     E_rot = B*Js*(Js+1) + (A-B)*K**2 - D*(Js*(Js+1))**2; # *hbar / hbar for solver

     # Note: since V and E_0_max ends up getting multiplied together,
     # 2 factors of epsilon_0 cancels
     # (since alpha is given in polarizability volume in Angstrom^3 =
     # 1/4pi epsilon_0 times polarizability )
     # Therefore, these factors are skipped.
     E_0_squared_max = 2*peak_intensity/c;
     
     if (not can_propagate_using_ODE):
         use_ODE = False;

     # If we don't propagate with the ODE method, instead
     # a matrix method is used. This method relies
     # on the diagonalized interaction matrix.
     if (not use_ODE or Jmax<=2):
        V0, V1, V2, eig, vec = interaction.interaction_matrix_div_E0_squared(Jmax,K,M,KMsign,delta_alpha,alpha_perp,diagonalize=True);
        # Divide by hbar for solver
        eig = eig*4*pi*1e-30/hbar
        psi_t = propagate(psi,t,E_rot,E_0_squared_max,sigma,eig,vec,xc_filename);
     else:
        V0, V1, V2 = interaction.interaction_matrix_div_E0_squared(Jmax,K,M,KMsign,delta_alpha,alpha_perp,diagonalize=False);
        V0 = V0*4*pi*1e-30/hbar; # Divide by hbar for solver
        V1 = V1*4*pi*1e-30/hbar;
        V2 = V2*4*pi*1e-30/hbar;
        psi_t = propagate_ODE(psi,t,E_rot,E_0_squared_max,sigma,V0,V1,V2,xc_filename,abstol,reltol);

     return psi_t;


def propagate(psi_0,time,E_rot,E_0_squared_max,sigma,eig,vec,xc_filename):
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
     if (len(E_rot) > 3):
        dt_max = min(2.4*sigma/150,1/(E_rot[-1]-E_rot[-3]));
     else:
        dt_max = min(2.4*sigma/150,1/E_rot[-1]);
     dt = time[1]-time[0];
     # If our given time step is larger than this, use a smaller time step
     if (dt > dt_max):
         scale = int(numpy.ceil(dt/dt_max));
         dt = dt/scale;
     else:
         scale = 1;

     customPulse=numpy.empty(time.size)
     use_custom_pulse=False


     if os.path.isfile(xc_filename):
         pulse_data = numpy.genfromtxt(xc_filename, delimiter=",");
         pulse_data[:,0] *= 1e-12;
         pulse_data[:,1] /= numpy.trapz(pulse_data[:,1],pulse_data[:,0]);
         fluence = E_0_squared_max/1e12; # 1e12 is because we normalize _after_ converting to seconds
         pulse_data[:,1] *= fluence;
         
         customPulse=linterp(pulse_data[:,0],pulse_data[:,1],time)
         use_custom_pulse=True
    #     scale=1

     expRot = numpy.exp(-1j*dt*E_rot);
     expRot2 = numpy.exp(-1j*dt*E_rot/2);

     psi_t = numpy.empty((len(time),len(psi_0)), dtype=numpy.complex)
     psi_t[0,:] = psi_0;
     vecT = numpy.ascontiguousarray(vec.T);

     if (libpropagation):
         #res = libpropagation.propagate_field(len(time), scale, len(psi_0), time[0], dt, E_0_squared_max, sigma, psi_t, alpha, t0, y0, eig, vec, vecT, expRot, expRot2);
         res = libpropagation.propagate_field(len(time),scale,len(psi_0),time[0],dt,E_0_squared_max,sigma,psi_t,eig,vec,vecT,expRot,expRot2,customPulse,use_custom_pulse);
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
                psi_t[i,:] = ((numpy.exp(-dt*E_0_squared*1j*eig))*vec).dot(vecT.dot(psi_t[i,:]));
                if (numpy.max(numpy.abs(psi_t[i,-2:])**2) > 1e-5):
                    raise RuntimeError("Basis size too small");
            psi_t[i,:] = expRot2*psi_t[i,:];
            i = i + 1;

     return psi_t;

def propagate_ODE(psi_0,time,E_rot,E_0_squared_max,sigma,V0,V1,V2,xc_filename,abstol=1e-8,reltol=1e-8):
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

def fieldfree_propagation(psi_0,t0,times,E_rot,Jmax,K,M,KMsign,D,do_cos2d=False):

     U, Udiag, U1, U2 = interaction.MeanCos2Matrix(Jmax,K,M,KMsign);
     cos2d = numpy.array([]);
     if (do_cos2d):
        U2d = U2dcalc.MeanCos2dMatrix(Jmax,K,M,KMsign);
     else:
        U2d = numpy.array([]);

     if (libpropagation and Jmax>0): # Call a C function instead of using numpy.
        psi = numpy.empty((len(times),Jmax+1),dtype=numpy.complex);
        cos2 = numpy.empty((len(times),),dtype=numpy.double);
        if (do_cos2d):
            cos2d = numpy.empty((len(times),),dtype=numpy.double);
    
        if (D > 0): #do we have centrifugal distortion?
            libpropagation.fieldfree_propagation(K,M,Jmax,psi_0,t0,len(times),times,E_rot,Udiag,U1,U2,psi,cos2,do_cos2d,U2d,cos2d, 1);
        else:
            libpropagation.fieldfree_propagation(K,M,Jmax,psi_0,t0,len(times),times,E_rot,Udiag,U1,U2,psi,cos2,do_cos2d,U2d,cos2d, 0);

        return psi,cos2,cos2d;

     phase = exp(-1j*numpy.outer(E_rot,(times-t0))); #should work with centrifugal dist
     psi = psi_0.reshape(Jmax+1,1)*phase;

     # We only need the diagonal in psi^H x U x psi.
     # since diag(AxB) = [sum A1j*Bj1, sum A2j*Bj2 ...] =
     # [sum A^Tj1*Bj1, sum A^Tj2*Bj2 ...] = A^T*B summed along the colums
     cos2 = numpy.real(numpy.conj(psi)*U.dot(psi)).sum(0);
     if (do_cos2d):
         cos2d = numpy.real(numpy.conj(psi)*U2d.dot(psi)).sum(0);

     psi=numpy.transpose(psi);

     return psi,cos2,cos2d; # psi[t][n] contains n'th component at time index t
