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
import argparse
import scipy
import dispatcher
import sys,os
import config
import utils

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate cos^2 trace of molecule in the |JKM> state.');

    parser.add_argument('molecule',type=str,help='Name of the molecule (e.g. CS2)');
    parser.add_argument('J',type=str,help='The J quantum number');
    parser.add_argument('K',type=str,help='The K quantum number');
    parser.add_argument('M',type=str,help='The M quantum number');
    parser.add_argument('moleculeconf',type=str,help='Config file for molecules');
    parser.add_argument('pulses',type=str,help='Laser pulse configuration file');
    parser.add_argument('--dt',nargs=1,type=float,default=[0.0],help='Temporal spacing in picoseconds (for plotting only)');
    parser.add_argument('-t',nargs=1,type=float,default=[-1],help='Propagate to this time (picoseconds)');
    parser.add_argument('--Jmax',nargs=1,default=[60],type=int,help='J of highest basis vector.');
    parser.add_argument('--Nshells',nargs=1,default=[1],type=int,help='Number of integration shells to use in focal volume averaging (1 for no focal volume averaging)');
    parser.add_argument('--probe_waist',nargs=1,default=[25.0],type=float,help='Probe waist in µm.');
    parser.add_argument('--cos2d',action="store_true",help='Calculate cos^2 theta 2d (requires the U2d helper program to be compiled, see Makefile)');
    parser.add_argument('--csv',action="store_true",help='Store results in a CSV file instead of a .npy file');
    parser.add_argument('-f','--filename',nargs=1,default=[""],type=str,help='Output filename. If using .npz format, .npz extension will be appended to the file name if it is not already there.');


    args = parser.parse_args();

    dt = args.dt[0]*1e-12;
    Nshells = args.Nshells[0];
    molecule = config.molecule(args.molecule,args.moleculeconf);
    pulses = config.load_laserpulses(args.pulses);
    J = int(args.J);
    K = int(args.K);
    M = int(args.M);
    Jmax = args.Jmax[0];
    if (dt <= 0):
        dt = utils.nice_time_step(molecule.B/(2*numpy.pi),Jmax)
    calculate_cos2d = args.cos2d;
    store_csv = args.csv;
    out_filename = args.filename[0];

    probe_waist = args.probe_waist[0]*1e-6;

    t_end = args.t[0]*1e-12;
    
    if (J < 0):
        raise RuntimeError('J must be positive');

   
    KMsign = numpy.sign(K*M);
    if (KMsign == 0):
        KMsign = 1;

    Kp = abs(K);
    Mp = abs(M);

    if (Mp > J or Kp > J):
        raise RuntimeError('K and M must be less than J');
    
    weight = 1;
    states = [(weight,J,Kp,Mp,KMsign)];
 
    do_psi_pulse = (Nshells == 1);

    t,cos2,cos2d,psi_pulse = dispatcher.dispatch(states,pulses,Jmax,Nshells,molecule,dt,t_end,probe_waist,calculate_cos2d,do_psi_pulse)
    psi_pulse = psi_pulse[0];

    if (do_psi_pulse):
        pdf = numpy.abs(psi_pulse)**2;
    
        Js = numpy.arange(0,Jmax+1);
        Jssq = Js**2;
        
        Javg = numpy.sum(Js*pdf,axis=1);
        Jsq_avg = numpy.sum(Jssq*pdf,axis=1);
    
        std = numpy.sqrt(Jsq_avg - Javg**2);
        cdf = numpy.cumsum(numpy.abs(psi_pulse)**2,axis=1);
        percentile_999 = numpy.argmax(cdf>=0.999,axis=1);


    attributes = [molecule.name,Jmax,J,K,M,dt,os.path.basename(args.pulses),Nshells,probe_waist];

    filename = out_filename;
    if (filename == ""):
        filename = ','.join([str(i) for i in attributes]) + ".npz";
        filename = "data/single_state_" + filename;

    try:
        os.mkdir("data");
    except OSError:
        pass;

    if (not store_csv):
     if (do_psi_pulse):
        numpy.savez(filename,t=t,cos2=cos2,cos2d=cos2d,Javg=Javg,std=std,percentile_999=percentile_999);
     else:
        numpy.savez(filename,t=t,cos2=cos2,cos2d=cos2d);
    else:
        if (out_filename == ""):
            filename = filename.replace("npz","csv");
        if (not do_psi_pulse):
            utils.save_to_csv(filename,t,cos2,cos2d);
        else:
            utils.save_to_csv(filename,t,cos2,cos2d,["<J>","std(J)","J_99.9%"],[Javg,std,percentile_999]);

    if (out_filename == ""):
        print("Saved trace in "+filename);



