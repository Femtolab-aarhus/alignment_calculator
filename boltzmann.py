#!/usr/bin/python2

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

import interaction
import numpy 
from numpy import pi,exp,sqrt,log,sin,cos
import sys


def thermal_JKM_ensemble(T,molecule,Jmax,percentile=0.999,anisotropy=1.0):
 
    # Returns a list of 5-tuples for each state in the ensemble.
    # Each tuple contains the weight, J, |K|, |M| and sign(KM).
    # The states with e.g. negative M give the same 
    # alignment trace as those with positive M (for K=0)
    # so only the positive M are chosen, and the weight is adjusted
    # accordingly.

    A = molecule.A;
    B = molecule.B;

    if (T == 0):
        if (molecule.even != 0):
            return [(1,0,0,0,1)];
        else:
            if (A != 0 and A-B < 0): # oblate top, K states energetically
                                     # favorable to M states.
                K = 1;
                w = 1/3.0;
                return [(w,1,K,0,1),(w,1,K,1,1),(w,1,K,1,-1)];
            else:
                K = 0;
                w = 1/3.0;
                return [(w,1,K,0,1), (2*w,1,K,1,1)];

    kb = 1.3806488e-23; # m^2kg/s^2 K
    hbar=1.05457173e-34;


    Js = numpy.arange(0,Jmax+1,1);

    if (molecule.A == 0): # Linear molecule
        Ks = numpy.array([0]);
    else:
        Ks = numpy.arange(0,Jmax+1,1);


    K, J = numpy.meshgrid(Ks,Js);
    E_rot = (B*J*(J+1) + (A-B)*K**2)*hbar;

    E_rot[K>J] = numpy.inf; # Remove unphysical K values by giving them infinite
                            # energy.
    
    # E_rot[J,K] now gives rotational energy of the |J,K> state

    g_m = (2*J+1); # M degeneracy

    z = exp(-E_rot/(kb*T));
    z[:,1:] *= 2; # for K > 0 -K,M gives same contribution to cos^2 as +K,-M,
    # Might want to throw some of these away because of nuclear spin etc.
    z[0::2,:] *= molecule.even;
    z[1::2,:] *= molecule.odd;
    Z = numpy.sum(g_m*z); # Partition function

    # Find the states that make up 99.9% of the partition function
    
    # Negative sign to sort in descending order
    order = (-g_m*z).flatten().argsort()

    z = z.flatten()[order];
    g_m = g_m.flatten()[order];
    Jf = J.flatten()[order];
    Kf = K.flatten()[order];

    part_z = numpy.cumsum(g_m*z);
    top = numpy.where(part_z/Z > percentile)[0][0]+1;

    #print(g_m*z/Z,file=sys.stderr);
    #print(Jf,file=sys.stderr);

    states = [];
    for weight,J,K in zip(z[:top],Jf[:top],Kf[:top]):
        # Don't calculate m<0 states, as these give
        # the same as m>0. Just give them double weight
        anisotropy_weight = anisotropy**numpy.arange(0,J+1);
        redist = numpy.ones_like(anisotropy_weight)*2;
        redist[0] = 1;
        anisotropy_weight /= numpy.sum(anisotropy_weight*redist)/(2*J+1);
        M = 0
        states.append([anisotropy_weight[M]*weight/Z,J,K,M,1]);
        for M in range(1,J+1):
            if (K > 0):
                states.append([anisotropy_weight[M]*weight/Z,J,K,M,1]);
                states.append([anisotropy_weight[M]*weight/Z,J,K,M,-1]);
            else:
                states.append([anisotropy_weight[M]*2*weight/Z,J,K,M,1]);

    w = [i[0] for i in states];
    total_weight = numpy.sum(numpy.array(w));
    if (total_weight > 1 or total_weight < percentile):
        raise RuntimeError("Something is wrong with the boltzmann weighting: percentile: {}, total weight: {}".format(percentile,total_weight));

    for state in states:
        state[0]/=total_weight; # Renormalize the weighting.

    return states;
