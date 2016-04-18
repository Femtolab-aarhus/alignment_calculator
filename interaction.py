#!/usr/bin/env python3

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
import random
import scipy
from scipy import linalg
from numpy import sqrt
import os,sys
import os.path
import tempfile
import utils


def A2(i): # Zare p63
     if (i <= 4):
          return 0;
     else:
          return 1/sqrt(i*(i-1.0)*(i-2.0)*(i-3.0)*(i-4.0));
   
def oddsign(i): #  i odd: -1, even: 1
     return -2*(i % 2 - 0.5);
     # Same as power(-1.0,i).


def interaction_matrix_div_E0_squared(Jmax,K,M,KMsign,delta_alpha,alpha_perp,diagonalize=True):
     ''' The interaction matrix <JKM|V|J'KM> divided by E_0^2
     (the electric field amplitude squared).
     E_0 is just a constant, so the same interaction matrix
     can be used for all laser pulses by just scaling this one
     with E_0^2 = 2*I/(c*epsilon_0) times a conversion factor
     between polarizability volume and polarizability.
     
     Input is the Jmax, the maximum J state,
     K and M for the state, along with KMsign, the relative sign between K and M,
     and delta_alpha and alpha_perp, the components of the polarizability
     tensor. 
     '''

     try:
         cache = interaction_matrix_div_E0_squared.cache;
     except AttributeError:
         cache = dict();
         interaction_matrix_div_E0_squared.cache = cache;
     try:
         return cache[(Jmax,K,M,KMsign,delta_alpha,alpha_perp,diagonalize)];
     except KeyError:
         pass;


     (U, U0, U1, U2) = MeanCos2Matrix(Jmax,K,M,KMsign);
     
     Jmin = max(K,M);

     tr = numpy.ones(Jmax+1); # Construct identity matrix,
     tr[0:Jmin] = 0;          # but with 0...Jmin zeroed

     V0 = -(delta_alpha*U0+alpha_perp*tr)/4;
     V1 = -delta_alpha*U1/4;
     V2 = -delta_alpha*U2/4;

     if (diagonalize):
         cat = numpy.concatenate;
         a = numpy.vstack((cat(([0,0],V2)),cat(([0],V1)),V0));
         eig,vec = linalg.eig_banded(a);
         vec = numpy.ascontiguousarray(vec);
         V = (V0, V1, V2, eig, vec);
     else:
         V = (V0, V1, V2);

     cache[(Jmax,K,M,KMsign,delta_alpha,alpha_perp,diagonalize)] = V;
     return V;


def MeanCos2Matrix(Jmax,K,M,KMsign):
     ''' Calculate <JKM|cos^2(theta)|J'KM>. '''
     
     Jmax = int(Jmax); # Avoids some overflow issues where e.g. Jmax is an int32
     K = int(K);
     M = int(M);
     KMsign = int(KMsign);

     try:
         cache = MeanCos2Matrix.cache;
     except AttributeError:
         cache = dict();
         MeanCos2Matrix.cache = cache;
     try:
         return cache[(Jmax,K,M,KMsign)];
     except KeyError:
         pass;


     # Wigner 3J symbols: CJJxM=(J J+x 2; M -M 0), see Zare p. 63
     
     U = numpy.zeros((Jmax+1,Jmax+1));
 
     Jmin = max(K,M);

     for J in range(Jmin, Jmax+1):
          CJJM = oddsign(J-M) * A2(2*J+3)*2*(3*M*M-J*(J+1));
          CJJK = oddsign(J-K) * A2(2*J+3)*2*(3*K*K-J*(J+1));
          
          U[J][J] = (2.0*(2*J+1)*oddsign(M-K)*CJJM*CJJK + 1.0)/3.0;

          if (J < Jmax):
               CJJ1M=oddsign(J+M+1)*A2(2*J+4)*2.0*M*sqrt(6*(J-M+1)*(J+M+1));
               CJJ1K=oddsign(J+K+1)*A2(2*J+4)*2.0*K*sqrt(6*(J-K+1)*(J+K+1));
               
               U[J][J+1] = 2.0*sqrt((2*J+1)*(2*J+3))*oddsign(M-K)*CJJ1M*CJJ1K*KMsign/3.0;
               U[J+1][J] = U[J][J+1];
          
          if (J < Jmax-1):
               CJJ2M=oddsign(J+M)*A2(2*J+5)*sqrt(6*(J+M+2)*(J+M+1)*(J-M+2)*(J-M+1));
               CJJ2K=oddsign(J+K)*A2(2*J+5)*sqrt(6*(J+K+2)*(J+K+1)*(J-K+2)*(J-K+1));
               U[J][J+2] = 2.0*sqrt((2*J+1)*(2*J+5))*oddsign(M-K)*CJJ2M*CJJ2K/3.0;
               U[J+2][J] = U[J][J+2];

     diag = numpy.ascontiguousarray(numpy.diag(U,0));
     U1 = numpy.ascontiguousarray(numpy.diag(U,1));
     U2 = numpy.ascontiguousarray(numpy.diag(U,2));
     cache[(Jmax,K,M,KMsign)] = (U, diag, U1, U2);
     return (U, diag, U1, U2);


def find_U2d_matrix_file(Jmax,K,M,KMsign):
    if (K != 0):
        raise NotImplementedError("K!=0 not implemented");

    beginning = '__U2d_array_';
    ending = '_'+str(K)+'_'+str(M)+ '_'+str(KMsign) + ".npy";
    names = os.listdir("./precomputed_U2d_matrices/");
    for name in names:
        if (name.startswith(beginning) and name.endswith(ending)):
            contents = name.split('_');
            Jmax_p = int(contents[4]);
            if (Jmax_p >= Jmax and name == beginning+str(Jmax_p)+ending):
                return name
    return False

def PrepareMeanCos2dMatrix(Jmax,K,M,KMsign):
    ''' Actual calculation of matrix elements. Call
        this function before calling MeanCos2dMatrix,
        as the latter will not do any actual calculation '''


    filename = find_U2d_matrix_file(Jmax,K,M,KMsign)
    if (not filename): # No suitable cos^2 theta 2d matrix found, calculate it
        if (utils.running_on_windows()):
            raise NotImplementedError("Generation of cos^2 theta 2D matrix elements currently not supported on the Windows platform.");

        signature = str(Jmax)+'_'+str(K)+'_'+str(M)+ '_'+str(KMsign);
        filename = './precomputed_U2d_matrices/__U2d_array_'+signature+'.npy';
    
        with tempfile.NamedTemporaryFile() as tmp:
            print("Generating U2d matrix Jmax,K,M,KMsign="+str((Jmax,K,M,KMsign))+".",file=sys.stderr);
            print("Saving result in: " + filename,file=sys.stderr);
            print("Temporary file: " + tmp.name,file=sys.stderr);
            command = './U2d '+tmp.name+' '+signature.replace("_", " ");
            if (os.system(command) != 0):
                raise RuntimeError("Failed to calculate U2d (command was "+command + ").");
            U2d = numpy.memmap(tmp.name,dtype=numpy.double,mode='r',shape=(Jmax+1,Jmax+1));
            numpy.save(filename,U2d);
    else:
        print("Reusing cos^2theta 2d matrix element file: "+filename,file=sys.stderr);


def MeanCos2dMatrix(Jmax,K,M,KMsign):
     ''' Calculate <JKM|cos^2(theta 2d)|J'KM>. '''
     
     try:
         cache = MeanCos2dMatrix.cache;
     except AttributeError:
         cache = dict();
         MeanCos2dMatrix.cache = cache;
     try:
         return cache[(Jmax,K,M,KMsign)];
     except KeyError:
         pass;

     filename = find_U2d_matrix_file(Jmax,K,M,KMsign)
     if (not filename):
         raise RuntimeError("Could not find any U2d matrix elements for Jmax, K, M, KMsign = " +str((Jmax,K,M,KMsign))+".");

     U2d = numpy.load("./precomputed_U2d_matrices/"+filename);
     U2d = numpy.ascontiguousarray(U2d[0:(Jmax+1),0:(Jmax+1)]);

     #for j in range(Jmax):
     #    for jp in range(j,Jmax,2):
     #        mult = 1.0;
     #        if (numpy.abs(U2d[j][jp]) < 1e-3):
     #            mult = 0.0;
     #        U2d[j][jp] *= mult;
     #        U2d[jp][j] *= mult;

     cache[(Jmax,K,M,KMsign)] = U2d;
     return U2d;


