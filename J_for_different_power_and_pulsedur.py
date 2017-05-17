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

# If you want to use this script, understand it, work on a copy of it
# and remove the following notice:

#print("You must edit this file before you use it!",file=sys.stderr)
# it is also a very good idea to run precalculate_matrix_elements.py before
# optimizing the temperature.
#import sys
#sys.exit(1)

import config
import propagation
import interaction
import numpy as np
from numpy import pi,exp,sqrt,log,sin,cos
import boltzmann
import dispatcher
import utils
import scipy.interpolate
import scipy.optimize


if __name__ == "__main__":
    # The program uses SI units internally.
    mol='OCS'

    p1=np.linspace(10,200,20)
    p2=np.linspace(250,500,6)
    p3=np.linspace(600,600,1)
    p=np.concatenate((p1,p2))
    p=np.concatenate((p,p3))
    print(p)
    Nshells = 1;
    # Replace I2 with the name of the molecule you run:
    molecule = config.molecule(mol,"./conf/molecules.conf");
    molecule.B=molecule.B/2.8
    waist=31e-6
    tau=12000e-15

    Jmax = 250;
    calculate_cos2d = False;
    dt = utils.nice_time_step(molecule.B/(2*numpy.pi),Jmax,for_cos2d=calculate_cos2d)

    probe_waist = 20e-6;
    percentile = 0.99;
    anisotropy = 1.0;

    t_end = 30e-12

    cache = dict();
    states = [(1,0,0,0,1)];
    extra_header = ["<J>","std(J)","J_99.9%","Probability coefficients"];
    delimiter = ","
    savepath='C:\Jmax_data/'
    for i in range(len(p)):
        INTENSITY=0.6*p[i]*1e-6/(tau*(waist**2))
        pulses = [config.laserpulse(INTENSITY,tau,0,waist=31e-6)];
        t,cos2,cos2d,psi = dispatcher.dispatch(states,pulses,Jmax,Nshells,molecule,dt,t_end,probe_waist,calculate_cos2d,do_psi_pulse=True)
        psi=psi[0]
        pdf = numpy.abs(psi)**2;  
        Js = numpy.arange(0,Jmax+1);
        Jssq = Js**2;
        Javg = numpy.sum(Js*pdf,axis=1);
        Jsq_avg = numpy.sum(Jssq*pdf,axis=1); 
        std = numpy.sqrt(Jsq_avg - Javg**2);
        cdf = numpy.cumsum(numpy.abs(psi)**2,axis=1);
        percentile_999 = numpy.argmax(cdf>=0.999,axis=1);
        psi_out=numpy.abs(psi[-1])**2
        extra_columns = [Javg,std,percentile_999,psi_out];
        filename=mol+'_p'+str(int(p[i]))+'_T'+str(tau*1e12)+'.csv'
        utils.save_to_csv(savepath+filename,t,cos2,cos2d,extra_header,extra_columns,delimiter);    





