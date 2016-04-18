#!/usr/bin/python3

import sys
import config
import propagation
import interaction
import numpy 
from numpy import pi,exp,sqrt,log,sin,cos
import boltzmann
import cos2_calc
import argparse
import time
import datetime
import clustercomputing
import ctypes
import tempfile

parser = argparse.ArgumentParser(description='Calculate transfer matrices for a molecule.');

parser.add_argument('molecule',type=str,help='Name of the molecule (e.g. CS2)');
parser.add_argument('moleculeconf',type=str,help='Config file for molecules');
parser.add_argument('pulses',type=str,help='Laser pulse configuration file');
parser.add_argument('authkey',type=str,help='Password for job server');
parser.add_argument('-T','--temperature',nargs=1,default=[1.],type=float,help='Rotational temperature in Kelvin');
parser.add_argument('--Jmax',nargs=1,default=[60],type=int,help='J of highest basis vector.');
parser.add_argument('--Nshells',nargs=1,default=[30],type=int,help='Number of integration shells to use in focal volume averaging (1 for no focal volume averaging)');
parser.add_argument('--probe_waist',nargs=1,default=[25.0],type=float,help='Probe waist in µm.');
parser.add_argument('--port',nargs=1,type=int,default=[2000],help='Listen on this port');

args = parser.parse_args();

Nshells = args.Nshells[0];
molecule = config.molecule(args.molecule,args.moleculeconf);
pulses = config.load_laserpulses(args.pulses);
Jmax = args.Jmax[0];
T = args.temperature[0];

probe_waist = args.probe_waist[0]*1e-6;

address = ('', args.port[0]);
authkey = bytes(args.authkey,'utf-8');

states = boltzmann.thermal_JKM_ensemble(T,molecule,Jmax,percentile=0.999);
focalvolume_weight,focal_pulses = cos2_calc.focal_volume_average(pulses,Nshells,probe_waist);
transfer_matrix_size = numpy.dtype(numpy.complex).itemsize*(Jmax+1)**2

manager,job_server = clustercomputing.start_server(address,authkey,timeout=10*60);

counter = 0;
save_args = dict();
ordered = dict();
for state in states:

    boltzmann_weight,J,K,M,KMsign = state;
    for n in range(Nshells):
        for pulse in focal_pulses[n]:
            Imax = pulse.I_max;
            FWHM = pulse.FWHM;
            t = 3*FWHM*numpy.array([-1,1]);
            if (propagation.load_transfer_matrix(molecule, Imax,FWHM,t[0],t[1],Jmax,K,M,KMsign) is None):
                job = [Imax,FWHM,t,Jmax,K,M,KMsign,molecule];
                order = (Imax,FWHM,t[0],t[-1],Jmax,K,M,KMsign,molecule);
                # Don't submit the same transfer matrix twice.
                try:
                    test = ordered[order];
                except KeyError:
                    counter += 1;
                    print(counter, round(counter*(Jmax+1)**2*8/1024/1024), "MB ", Jmax, FWHM, K, M);
                    job_id = job_server.submit_job(job);
            
                    save_args[job_id] = [molecule,Imax,FWHM,t[0],t[1],Jmax,K,M,KMsign];
                    ordered[order] = True;



job_server.close();
print("Submitted "+str(counter)+" jobs. Direct workers to port " + str(address[1]));
for c in range(counter):
    job_id,transfer_matrix = job_server.get_result();
    
    args = save_args[job_id];
    args.append(transfer_matrix);
    I_max,FWHM,t0,t1,Jmax,K,M,KMsign,mat= args[1:]

    print("Saving Jmax,K,M,KMs,I_max,FWHM = {:s},{:s},{:s},{:s},{:s},{:s}".format(*[str(i) for i in [Jmax,K,M,KMsign,I_max,FWHM]]));
    print(str(c/counter*100) + " % done.");
    propagation.save_transfer_matrix(*args);
    
    del save_args[job_id];

manager.shutdown();




