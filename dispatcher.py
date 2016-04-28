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

import kill_library_multithreading # must be done first

import numpy 
from numpy import pi
import time
import datetime
import trace_backend
import multiprocessing
import sys,os


def dispatch(states,pulses,Jmax,Nshells,molecule,dt,t_end,probe_waist,calculate_cos2d,do_psi_pulse=False):
       
    if (t_end < 0):
        B = molecule.B;
        revival = 1/(2*B/(2*pi));
        t_end = 1.1*revival;

    cos2s = [];
    cos2ds = [];
    psis = [];

    focalvolume_weight,focal_pulses = trace_backend.focal_volume_average(pulses,Nshells,probe_waist);
 
    dispatcher_started_time = time.time();
    with multiprocessing.Pool() as p:

        asyncs = [];
        for state in states:
            boltzmann_weight,J,K,M,KMsign = state;
            for n in range(Nshells):
                asyncs.append(p.apply_async(trace_backend.cos2_trace,(J,K,M,KMsign,Jmax,molecule,focal_pulses[n],dt,t_end,calculate_cos2d,do_psi_pulse)));

        num = 0;
        for state in states:
            boltzmann_weight,J,K,M,KMsign = state;
            print("Waiting for cos2 for J,K,M,KMsign = {:d},{:d},{:d},{:d}.".format(J,K,M,KMsign),file=sys.stderr);
            cos2 = 0;
            cos2d = 0;
            begin = time.time()
            for n in range(Nshells):
                t,cos2_n,cos2d_n,psi_n = asyncs.pop(0).get();
                cos2 += focalvolume_weight[n]*cos2_n;
                cos2d += focalvolume_weight[n]*cos2d_n;
            cos2s.append(boltzmann_weight*cos2);
            cos2ds.append(boltzmann_weight*cos2d);
            psis.append(psi_n);
            num = num + 1;
            timedelta = timedelta_round(datetime.timedelta(seconds=time.time()-begin),1);
            print(("Time delta: "+timedelta).ljust(40) + "{:.1f}% done.".format(num/len(states)*100).rjust(40)  ,file=sys.stderr);

    print("",file=sys.stderr);
    timedelta = datetime.timedelta(seconds=time.time()-dispatcher_started_time);
    timedelta = timedelta_round(timedelta,1);
    print("Job completed in total time: "+timedelta,file=sys.stderr);

    cos2 = numpy.array(cos2s);
    cos2 = numpy.sum(cos2,0);
    cos2d = numpy.array(cos2ds);
    cos2d = numpy.sum(cos2d,0);

    return t,cos2,cos2d,psis

def timedelta_round(td,digits):
    td = str(td).split(".");
    if (len(td) > 1):
        fraction = "{:.0f}".format(round(float(td[1])/1e6,digits)*10**digits)
        td = td[0]+"."+fraction;
    else:
        td = td[0];
    return td;

