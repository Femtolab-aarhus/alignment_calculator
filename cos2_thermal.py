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

import sys,os
import config
import propagation
import interaction
import numpy 
from numpy import pi,exp,sqrt,log,sin,cos
import boltzmann
import argparse
import time
import datetime
import multiprocessing
import dispatcher
import utils


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate cos^2 trace of molecule.');

    parser.add_argument('molecule',type=str,help='Name of the molecule (e.g. CS2)');
    parser.add_argument('moleculeconf',type=str,help='Config file for molecules');
    parser.add_argument('pulses',type=str,help='Laser pulse configuration file');
    parser.add_argument('--dt',nargs=1,type=float,default=[0.0],help='Temporal spacing in picoseconds (for plotting only)');
    parser.add_argument('-t',nargs=1,type=float,default=[-1],help='Propagate to this time (picoseconds)');
    parser.add_argument('-T','--temperature',nargs=1,default=[1.],type=float,help='Rotational temperature in Kelvin');
    parser.add_argument('--Jmax',nargs=1,default=[60],type=int,help='J of highest basis vector.');
    parser.add_argument('--Nshells',nargs=1,default=[30],type=int,help='Number of integration shells to use in focal volume averaging (1 for no focal volume averaging)');
    parser.add_argument('--probe_waist',nargs=1,default=[25.0],type=float,help='Probe waist in µm.');
    parser.add_argument('--percentile',nargs=1,default=[0.999],type=float,help='Resolve this percentile of the Boltzmann distribution. Must be between 0 and 1 exclusive.');
    parser.add_argument('--anisotropy',nargs=1,default=[1.0],type=float,help='An Even Lavie vale tends to skew the equipartitioned initial state towards either high or low M. Multiply the Boltzmann weight by this factor raised to the power M. Values lower than 1 means low M becomes more propable. Values higher than 1 means high M becomes more probable.');
    parser.add_argument('--cos2d',action="store_true",help='Calculate cos^2 theta 2d (requires the U2d helper program to be compiled, see Makefile)');
    parser.add_argument('--csv',action="store_true",help='Store results in a CSV file instead of a .npy file');
    parser.add_argument('-f','--filename',nargs=1,default=[""],type=str,help='Output filename. If using .npz format, .npz extension will be appended to the file name if it is not already there.');


    args = parser.parse_args();

    dt = args.dt[0]*1e-12;
    Nshells = args.Nshells[0];
    molecule = config.molecule(args.molecule,args.moleculeconf);
    pulses = config.load_laserpulses(args.pulses);
    Jmax = args.Jmax[0];
    if (dt <= 0):
        dt = utils.nice_time_step(molecule.B/(2*numpy.pi),Jmax)
    T = args.temperature[0];
    calculate_cos2d = args.cos2d;
    store_csv = args.csv;
    out_filename = args.filename[0];

    probe_waist = args.probe_waist[0]*1e-6;
    percentile = args.percentile[0];
    if (percentile >= 1 or percentile <=0):
        raise ValueError("Percentile must be between 0 and 1, both not included.");
    anisotropy = args.anisotropy[0];
    if (anisotropy < 0):
        raise ValueError("Anisotropy factor must be nonnegative.");

    t_end = args.t[0]*1e-12;

    states = boltzmann.thermal_JKM_ensemble(T,molecule,Jmax,percentile,anisotropy);

    t,cos2,cos2d,psi = dispatcher.dispatch(states,pulses,Jmax,Nshells,molecule,dt,t_end,probe_waist,calculate_cos2d,do_psi_pulse=False)

    attributes = [molecule.name,Jmax,dt,os.path.basename(args.pulses),Nshells,T,probe_waist];

    filename = out_filename;
    if (filename == ""):
        filename = "data/"+','.join([str(i) for i in attributes]) + ".npz";
    
    meta = dict();
    meta['molecule'] = molecule;
    meta['Jmax'] = Jmax;
    meta['dt'] = dt;
    meta['pulses'] = pulses;
    meta['Nshells'] = Nshells;
    meta['temperature'] = T;
    meta['probe_waist'] = probe_waist;

    try:
        os.mkdir("data");
    except OSError:
        pass;

    if (not store_csv):
        numpy.savez(filename,t=t,cos2=cos2,cos2d=cos2d,meta=meta);
    else:
        if (out_filename == ""):
            filename = filename.replace("npz","csv");
        utils.save_to_csv(filename,t,cos2,cos2d);

    if (out_filename == ""):
        print("Saved trace in "+filename);


