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
from numpy import sqrt,log,pi
import itertools
import argparse

parser = argparse.ArgumentParser(description='');
parser.add_argument('intensity',type=str,help='Intensity in mW (e.g. 100, or 100:50:200');
parser.add_argument('FWHM',type=str,help='FWHM in fs.');

args = parser.parse_args();
I = args.intensity.split(":");
FWHM = args.FWHM.split(":");
I = [float(i) for i in I];
FWHM = [float(i) for i in FWHM];
if (len(I) == 1):
    I = numpy.array(I)*1e-3;
else:
    I = numpy.arange(I[0],I[2],I[1])*1e-3;

if (len(FWHM) == 1):
    FWHM = numpy.array(FWHM)*1e-15;
else:
    FWHM = numpy.arange(FWHM[0],FWHM[2],FWHM[1])*1e-15;

# Intensity as measured in watt in the lab
#I = numpy.arange(100,160,10)*1e-3;
#I = numpy.array([100])*1e-3;
# Pulse duration (FWHM) as measured in the lab [seconds]
#FWHM = numpy.arange(100,400,50)*1e-15;
#FWHM = numpy.array([200])*1e-15;
# Kick waist as measured in the lab (does not depend on I or FWHM(?))
kick_waist = 35e-6; #meter

repetition_rate = 1e3; # Hz

pulse_energy = I/repetition_rate;

E_grid,Tau_grid = numpy.meshgrid(pulse_energy,FWHM);

peak_intensity = 4*sqrt(log(2))/pi**(3.0/2)*E_grid/(kick_waist**2*Tau_grid);
peak_intensity = numpy.around(peak_intensity,-16);

I_max = peak_intensity.flatten();
FWHM = Tau_grid.flatten();


print("#Units: (SI)");
print("# I_max: W/m^2");
print("# FWHM, t: seconds");
print("# waist: meters");
print("[pulses]");
print("I_max = " + ', '.join([str(I) for I in I_max]).replace('+',''));
print("FWHM = " + ', '.join([str(T) for T in FWHM]));
print("t = " + ','.join(["0" for i in I_max]));
print("waist = " + ','.join([str(kick_waist) for i in I_max]));




