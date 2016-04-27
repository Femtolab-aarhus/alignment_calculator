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


# Example on how to fit 2d alignment traces via the temperature.
# It can be done faster without interpolation by only propagating directly
# to the data points. That's just a bit more cumbersome..

import sys
print("You must edit this file before you use it!",file=sys.stderr)
# it is also a very good idea to run precalculate_matrix_elements.py before
# optimizing the temperature.
sys.exit(1)

import config
import propagation
import interaction
import numpy 
from numpy import pi,exp,sqrt,log,sin,cos
import boltzmann
import dispatcher
import utils
import scipy.interpolate
import scipy.optimize


if __name__ == "__main__":
    # The program uses SI units internally.

    data = numpy.genfromtxt("./your_data.csv");
    t_data = data[:,0]*1e-12; # ps to seconds
    cos2d_data = data[:,1];
    filt = (t_data >=-1); # Remove data that the program does not calculate.
    t_data = t_data[filt];
    cos2d_data = cos2d_data[filt];

    Nshells = 20;
    # Replace I2 with the name of the molecule you run:
    molecule = config.molecule("I2","./conf/molecules.conf");
    pulses = [config.laserpulse(INTENSITY,FWHM_DURATION,0,waist=35e-6)];
    Jmax = 40;
    calculate_cos2d = True;
    dt = utils.nice_time_step(molecule.B/(2*numpy.pi),Jmax,for_cos2d=calculate_cos2d)

    probe_waist = 20e-6;
    percentile = 0.999;
    anisotropy = 1.0;

    t_end = t_data[-1];

    cache = dict();

    def sum_sq_err(T):
        print(T)
        if (T < 0):
            return numpy.inf;

        states = boltzmann.thermal_JKM_ensemble(T,molecule,Jmax,percentile,anisotropy);
        t,cos2,cos2d,psi = dispatcher.dispatch(states,pulses,Jmax,Nshells,molecule,dt,t_end,probe_waist,calculate_cos2d,do_psi_pulse=False)
        
        f = scipy.interpolate.interp1d(t,cos2d,kind='linear');
        cos2d_calc = f(t_data);
        error = numpy.sum((cos2d_calc-cos2d_data)**2);
        return error;
    
    T_guess = numpy.array([1]); # kelvin
    T_fit = scipy.optimize.fmin(sum_sq_err,T_guess);
    print(T_fit);



