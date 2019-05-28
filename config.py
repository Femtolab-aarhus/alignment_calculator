
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

from numpy import pi
import configparser
import sys
from functools import total_ordering

class molecule(object):
     
     def __init__(s, name, filename=None):
          s.name = name;
          if (filename is not None):
               s.load_molecule(filename);
     
     def load_molecule(s,filename):
         
          s.filename=filename;

          parser = configparser.ConfigParser();

          parser.read(filename);
        
          s.A_raw = parser.getfloat(s.name,"A");
          s.B_raw = parser.getfloat(s.name,"B");
          s.D_raw = parser.getfloat(s.name,"D");

          s.A = 2*pi*s.A_raw*1e9; # Reexpress as
          s.B = 2*pi*s.B_raw*1e9; # Hz angular frequency
          s.D = 2*pi*s.D_raw*1e3; #D is in kHz 
                # E = hbar*B*j*(j+1) -> the time evolution exp(-iE/hbar*dt) 
                # becomes exp(-iB*j*(j+1)*dt)
          s.alpha_par_volume = parser.getfloat(s.name,"alpha_par_volume");
          s.alpha_perp_volume = parser.getfloat(s.name,"alpha_perp_volume");
          s.delta_alpha = s.alpha_par_volume - s.alpha_perp_volume;
          s.even = parser.getfloat(s.name,"even");
          s.odd = parser.getfloat(s.name,"odd");
     
          #s.hbar = parser.getfloat(name,"hbar");
          #s.speed_of_light = parser.getfloat(name,"speed_of_light");

     def __repr__(s):
          return "molecule(\"" + str(s.name) + "\", \""+str(s.filename)+"\")";

     def __str__(s):
          return s.__repr__();



@total_ordering
class laserpulse(object):
    def __init__(s,I_max,FWHM,t,waist=None):

            s.I_max = float(I_max);
            s.FWHM = float(FWHM);
            s.t = float(t);
            if (waist):
                s.waist = float(waist);
    def __eq__(s,other):
        return s.t == other.t;
    def __lt__(s,other):
        return s.t<other.t;

    def __repr__(s):
        try:
            return "laserpulse({:s}, {:s}, {:s}, {:s})".format(str(s.I_max),str(s.FWHM),str(s.t),str(s.waist));
        except AttributeError:
            return "laserpulse({:s}, {:s}, {:s})".format(str(s.I_max),str(s.FWHM),str(s.t));

    def __str__(s):
        return s.__repr__();


def load_laserpulses(filename):
    parser = configparser.ConfigParser();

    parser.read(filename);

    I_max = parser.get("pulses","I_max").split(",");
    FWHM = parser.get("pulses","FWHM").split(",");
    t = parser.get("pulses","t").split(",");
    try:
        waist = parser.get("pulses","waist").split(",");
    except configparser.NoOptionError:
        waist = len(I_max)*[None];
    
    if (len(I_max) != len(FWHM) or len(FWHM) != len(t) or len(I_max)!=len(waist)):
        print("Malformed pulse file",file=sys.stderr);
        sys.exit(1);
    
    pulses = [];
    for i in range(len(I_max)):
        pulses.append(laserpulse(I_max[i],FWHM[i],t[i],waist[i]));

    return pulses;

