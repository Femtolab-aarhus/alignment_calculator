#!/usr/bin/python3 

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
import utils
import ctypes
import ctypes as ct
from numpy.ctypeslib import ndpointer
import scipy.integrate
import scipy.special
import sys,os
from itertools import product

matrix_element_file = "./precomputed_matrix_elements.npy";
U2dlib = None;
try:
    if (utils.running_on_windows()):
        U2dlib = ctypes.CDLL("./libU2d.dll");
    else:
        U2dlib = ctypes.CDLL("./libU2d.so");
    real_ndptr = ndpointer(flags=("CONTIGUOUS",),dtype=numpy.double);

    try: # May not be implemented.
        U2dlib.expand_U2d.restype = None;
        U2dlib.expand_U2d.argtypes = (ct.c_int,real_ndptr);
    except:
        pass;

    U2dlib.populate_U2d.restype = None;
    U2dlib.populate_U2d.argtypes = (ct.c_int,ct.c_int,ct.c_int,ct.c_int,real_ndptr,real_ndptr);
    
    # Disable multithreading in U2dlib. We are allready parallelized.
    try: 
        U2dlib.omp_set_num_threads(1); # May fail if library compiled witout omp
    except:
        pass;

except OSError:
    print("Not using C matrix element library. Compile it with make.",file=sys.stderr);

def set_num_threads(num):
    try:
        if (num>0):
            U2dlib.omp_set_num_threads(num);
    except:
        pass;

def MeanCos2dMatrix(Jmax,K,M,KMsign):
    ''' Calculate <JKM|cos^2(theta 2d)|J'KM>. '''
    
    if (K < 0 or M < 0):
        raise ValueError("K and M must both be nonnegative.");

    try:
        cache = MeanCos2dMatrix.cache;
    except AttributeError:
        cache = dict();
        MeanCos2dMatrix.cache = cache;
    try:
        return cache[(Jmax,K,M,KMsign)];
    except KeyError:
        pass;

    # See if we have some precalculated matrix elements available..
    U2d = filesystem_cache_lookup(Jmax,K,M,KMsign);

    if (U2d is None):
        U2d = actually_calculate_U2d(Jmax,K,M,KMsign);
    
    cache[(Jmax,K,M,KMsign)] = U2d;
    return U2d;

 
def filesystem_cache_lookup(Jmax,K,M,KMsign):
    ''' See if there is already a precomputed cache available before runtime '''

    try:
        cache = filesystem_cache_lookup.cache
    except AttributeError:
        cache = dict();
        cache["FileNotFound"] = False;
        filesystem_cache_lookup.cache = cache;
    
    if (cache["FileNotFound"]):
        return None;

    try:
        cache_file = cache["array"];
    except KeyError:
        try:
            # Use memory mapping to avoid loading the entire file into memory.
            # Chances are, we're only going to need a thin slice of it.
            cache_file = numpy.load(matrix_element_file,mmap_mode='r');
            cache["array"] = cache_file;
        except FileNotFoundError:
            cache["FileNotFound"] = True; 
            return None

    try:
        (kms, cKmax, cMmax, cJmax, cJmaxp) = cache_file.shape
        cMmax -= 1;
        cKmax -= 1;
        cJmax -= 1;
    except ValueError:
        print("precomputed_matrix_elements malformed.",file=sys.stderr);
        cache["FileNotFound"] = True;
        return None

    if (M <= cMmax and K <= cKmax and Jmax <= cJmax):
        if (KMsign == 1):
            KMsign = 0;
        KMsign *=-1;
        return numpy.ascontiguousarray(cache_file[KMsign,K,M,0:(Jmax+1),0:(Jmax+1)]);
    else:
        return None

def cache_size():
    try:
        cache_file = numpy.load(matrix_element_file,mmap_mode='r');
    except FileNotFoundError:
        return None;

    try:
        (kmsgn, Kmax, Mmax, Jmax, Jmaxp) = cache_file.shape
        Kmax -= 1;
        Mmax -= 1;
        Jmax -= 1;
    
        return (Jmax,Kmax,Mmax);
    except:
        return None;

def drop_cache():
    try:
        os.unlink(matrix_element_file);
    except:
        pass;


def actually_calculate_U2d(Jmax,K,M,KMsign):

    if (U2dlib is None):
        raise FileNotFoundError("libU2d not found. This is required for calculating cos^2 theta 2d matrix elements.");

    try:
        expansion_coefficients = actually_calculate_U2d.expansion_coefficients;
    except AttributeError:
        # Expand cos^2 theta 2d in legendre polynomials P_l(u)
        # for l = 0,1,2,...,lmax
        lmax = 100; # min(100,2*Jmax); logic error if first call jmax=2, second jmax=200
        expansion_coefficients = numpy.empty(lmax+1);
        try:
            U2dlib.expand_U2d(lmax,expansion_coefficients);
        except AttributeError:
            # that required Gnu Scientific Library. If that was not
            # found, calculate them slowly using python
            expansion_coefficients_in_python(lmax,expansion_coefficients);
        # Store the coefficients for later.
        actually_calculate_U2d.expansion_coefficients = expansion_coefficients;

    U2d = numpy.empty((Jmax+1,Jmax+1));
    U2dlib.populate_U2d(Jmax,KMsign*K,M,len(expansion_coefficients)-1,expansion_coefficients,U2d);
    
    return U2d;

def expansion_coefficients_in_python(lmax,a):
    # Calculate cos^2 theta 2d expansion coefficients in python instead
    # of using the C library.
    
    print("C call not available. Calculating expansion coefficients in Python instead.\nThis is slow and can take a while...",file=sys.stderr);

    tol = 1e-13;

    def inner(phi,u):
        if (u == 0):
            return 0;
        usq = u*u;
        return usq/(usq+(1-usq)*numpy.sin(phi)**2);
    
    def outer(u,l):
        theta = numpy.arccos(u);
        res = scipy.integrate.quad(inner,0,2*numpy.pi,args=(u,),epsabs=tol,epsrel=tol,limit=409600)[0];
        return numpy.real(scipy.special.sph_harm(0,l,0,theta))*res;

    a[:] = 0;
    for l in range(0,lmax+1,2):
        a[l] = scipy.integrate.quad(outer,-1.0,1.0,args=(l,),epsabs=tol,epsrel=tol,limit=409600)[0];

    print("Done expanding cos^2 theta 2d.",file=sys.stderr);

def precalculate_matrix_elements(Jmax,Kmax,Mmax):

    a = numpy.zeros((2,Kmax+1,Mmax+1,Jmax+1,Jmax+1));

    KMsigns = (1,-1);
    Ks = range(0,Kmax+1);
    Ms = range(0,Mmax+1);
     
    for (KMsign,K,M) in product(KMsigns,Ks,Ms):
        if (KMsign == -1 and (K == 0 or M == 0)):
            continue;
        sign_idx = KMsign;
        if (sign_idx==1):
            sign_idx = 0;
        else:
            sign_idx=1;
        a[sign_idx,K,M,:,:] = MeanCos2dMatrix(Jmax,K,M,KMsign)[:,:];

    print("Saving matrix elements to: ",matrix_element_file);
    numpy.save(matrix_element_file,a);

