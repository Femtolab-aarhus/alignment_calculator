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
from numpy import pi,exp
import propagation
import interaction
import copy
import scipy
import dispatcher
import sys,os
import config
import utils
import time
import U2dcalc

def cos2_trace(J,K,M,KMsign,Jmax,molecule,laserpulses,dt,t_end,do_cos2d,do_psi_pulse=False):
    laserpulses.sort();
    window = 3/2.0;

    psi = numpy.zeros(Jmax+1,dtype=numpy.complex);
    psi[J] = 1;

    times = [];
    cos2 = [];
    cos2d = [];
    psis = [];


    U, U0, U1, U2 = interaction.MeanCos2Matrix(Jmax,K,M,KMsign)
    if (do_cos2d):
        U2d = U2dcalc.MeanCos2dMatrix(Jmax,K,M,KMsign);

    last_t = laserpulses[0].t-window*laserpulses[0].FWHM-max(1e-12,5*dt); # Start 1 ps before first pulse

    A = molecule.A;
    B = molecule.B;
    Js = numpy.arange(0,Jmax+1,1);
    E_rot = (B*Js*(Js+1) + (A-B)*K**2); # * hbar but / hbar for input to fieldfree
    for pulse in laserpulses:

        t = pulse.t;
        FWHM = pulse.FWHM;
        I_max = pulse.I_max;

        if (last_t > t-window*FWHM):
            raise RuntimeError("Pulses are not well enough separated.");

        num_steps = max(2,numpy.ceil((t-window*FWHM-last_t)/dt))
        times_before = numpy.linspace(last_t,t-window*FWHM,num_steps);
        # Propagate between pulses
        psi_before,cos2_before,cos2d_before = propagation.fieldfree_propagation(psi,last_t,times_before,E_rot,Jmax,K,M,KMsign,do_cos2d);

        num_steps = max(2,numpy.ceil(2*window*FWHM/dt));
        if (do_psi_pulse):
            num_steps = max(20,num_steps);

        integration_time = numpy.linspace(-window*FWHM,window*FWHM,num_steps);
        transfer = propagation.transfer_KM(psi_before[-1],K,M,KMsign,Jmax,I_max,FWHM,integration_time,molecule);
        S = sum(numpy.abs(transfer[-1,:])**2)
        if (numpy.abs(1-S) > 0.001):
            raise RuntimeError("Norm not preserved! "+ str((J,K,M,S)));

        times.append(times_before[1:-1]);
        cos2.append(cos2_before[1:-1]);
        cos2d.append(cos2d_before[1:-1]);

        times.append(integration_time+t)
        cos2_pulse = numpy.empty((len(integration_time),))
        cos2d_pulse = numpy.empty((len(integration_time),))
        for i in range(len(integration_time)):
            psi = transfer[i,:];
            cos2_pulse[i] = numpy.real(numpy.conj(psi).dot(U.dot(psi)))
            if (do_cos2d):
                cos2d_pulse[i] = numpy.real(numpy.conj(psi).dot(U2d.dot(psi)));
        cos2.append(cos2_pulse)
        if (do_cos2d):
            cos2d.append(cos2d_pulse);

        if (do_psi_pulse):
            psis.append(psi_before[1:-1,:]);
            psis.append(transfer);

        psi = transfer[-1,:];
        last_t = t + window*FWHM;

    # End looping over laser pulses. Now propagate to t_end.

    times_after = numpy.arange(last_t,last_t+t_end,dt);
    psi_after,cos2_after,cos2d_after = propagation.fieldfree_propagation(psi,last_t,times_after,E_rot,Jmax,K,M,KMsign,do_cos2d);
    times.append(times_after[1:]);
    cos2.append(cos2_after[1:]);
    if (do_cos2d):
        cos2d.append(cos2d_after[1:]);

    times = numpy.concatenate(times);
    cos2 = numpy.concatenate(cos2);
    cos2d = numpy.concatenate(cos2d);

    if (do_psi_pulse):
        psis = numpy.concatenate(psis);
    return times,cos2,cos2d,psis


def focal_volume_average(laserpulses,Nshells,probe_waist):
    
    pulses = [];

    laserpulses.sort();

    rmax = 1.7*probe_waist;
    dr = rmax/Nshells;

    n = numpy.arange(0,Nshells);
    r = (2.0*n+1.0)*dr/2.0;
    Volume = 2*r*dr;
    Pdiss = exp(-2.0*r**2/(probe_waist**2));
    Weight = Volume*Pdiss;

    Weight[0] = dr*dr;
    Weight /= numpy.sum(Weight);

    shell_pulses = copy.deepcopy(laserpulses);
    pulses.append(shell_pulses);

    for n in range(1,Nshells):
        shell_pulses = copy.deepcopy(laserpulses);
        for i in range(len(laserpulses)):
            waist = laserpulses[i].waist;
            shell_pulses[i].I_max = laserpulses[i].I_max*exp(-2.0*r[n]**2/(waist**2));
        pulses.append(shell_pulses);

    return Weight,pulses;

