#!/usr/bin/python3

import argparse
import sqlite3
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy
import propagation
import config

parser = argparse.ArgumentParser(description='Plot real and imaginary part of transfer matrix along with abs^2');
parser.add_argument('molecule',type=str,help='Name of the molecule (e.g. CS2)');
parser.add_argument('molecule_conf',type=str,help='Name of the molecules config file ');
parser.add_argument('I_max',type=float,help='Peak intensity (TW/cm^2)');
parser.add_argument('Jmax',type=int,help='Maximum J');
parser.add_argument('K',type=int,help='K');
parser.add_argument('M',type=int,help='M');
parser.add_argument('KMsign',type=int,help='KMsign');
parser.add_argument('FWHM',type=float,help='FWHM (in ps) of laser pulse');
parser.add_argument('window',type=float,help='Number of FWHM on either side of pulse');

args = parser.parse_args();

I_max = args.I_max;
Jmax = args.Jmax;
K = args.K;
M = args.M;
KMsign = args.KMsign;
FWHM = args.FWHM;
window = args.window;

molecule = config.molecule(args.molecule,args.molecule_conf);

t = window*FWHM*numpy.array([-1,1]);
transfer_matrix = propagation.get_transfer_matrix(I_max,FWHM,t,Jmax,K,M,KMsign,molecule);
if (transfer_matrix is None):
    print("Transfer matrix not stored. Calculating it.");
    transfer_matrix = propagation.compute_transfer_matrix(I_max,FWHM,t,Jmax,K,M,KMsign,molecule);
    propagation.save_transfer_matrix(molecule,I_max,FWHM,t[0],t[-1],Jmax,K,M,KMsign,transfer_matrix);

cmap = cm.copper_r
#cmap = cm.prism_r
cmap = cm.ocean_r


plt.matshow(numpy.real(transfer_matrix),cmap=cmap);
plt.title('Real part');
plt.colorbar();
plt.matshow(numpy.imag(transfer_matrix),cmap=cmap);
plt.title('Imaginary part');
plt.colorbar();
plt.matshow(numpy.abs(transfer_matrix)**2,cmap=cmap);
plt.title('Norm squared');
plt.colorbar();
plt.matshow(numpy.angle(transfer_matrix,True),cmap=cmap);
plt.title('Angle');
plt.colorbar();



plt.show();
