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
import platform


def save_to_csv(filename,t,cos2,cos2d,extra_headers=[],extra_columns=[],delimiter=" "):
    csvfile = numpy.empty((len(t),3));
    csvfile[:,0] = t/1e-12;
    csvfile[:,1] = cos2;
    if (len(cos2d)==len(t)):
        csvfile[:,2] = cos2d;
    else:
         csvfile[:,2] = numpy.nan;

    with open(filename,"wb") as f:
        header = "Time[ps], cos^2 theta, cos^2 theta 2d";
        for extra in extra_headers:
            header += " "+extra;
        if (len(extra_columns) > 0):
            extra = numpy.empty((csvfile.shape[0],len(extra_columns)))*numpy.nan;
            for i in range(len(extra_columns)):
                extra[:len(extra_columns[i]),i] = extra_columns[i]
            csvfile = numpy.concatenate((csvfile,extra),axis=1)
        numpy.savetxt(f,csvfile,delimiter=delimiter,header=header);


def running_on_windows():
    return platform.system().lower().startswith('win');

def nice_time_step(B,Jmax):
    ''' Chooses a nice time step based on the molecule and basis size '''
    # The highest frequency is the beat between
    # J = 0 and J = Jmax. However, this is never realized because
    # of the Delta J = 2 selection rule for the pulse.
    # The highest frequency components of any observable
    # is therefore the frequency of the beat between Jmax and Jmax-2
    # 
    # Although K>0 shifts the energy, the pulse does not change K,
    # so this offset is the same for all states.
    # Therefore, this calculation does not depend on A.
    # 
    samples_per_period = 10; # 2 for Nyquist criterion
    # Note: B given in Hz. The oscillation frequency
    # is simply the energy.
    E_max = B*Jmax*(Jmax+1);
    E_second = B*(Jmax-2)*(Jmax-1);
    
    frequency = E_max-E_second;
    sampling_frequency = frequency*samples_per_period;
    dt = 1.0/(sampling_frequency);
    
    return dt;


