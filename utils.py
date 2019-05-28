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
import numpy as np


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
def running_on_mac():
    return platform.system().lower().startswith('darwin');

def library_file_ending():
    if (running_on_windows()):
        return "dll"
    if (running_on_mac()):
        return "dylib"
    else:
        return "so"

def nice_time_step(B,Jmax,for_cos2d=False):
    ''' Chooses a nice time step based on the molecule and basis size '''
    # The highest frequency is the beat between
    # J = 0 and J = Jmax. However, this is never realized because
    # of the Delta J = 2 selection rule for the observable.
    # The highest frequency components of cos^2 theta 2d
    # is therefore the frequency of the beat between Jmax and Jmax-2
    #
    # For cos^2 theta_2d, there are nonzero components out to at least J=100.
    # So the rule here is beats between J and J±100 at least.
    # However, these are very small and J±4 should suffice.
    #
    # Although K>0 shifts the energy, the pulse does not change K,
    # so this offset is the same for all states.
    # Therefore, this calculation does not depend on A.
    #
    samples_per_period = 10; # 2 for Nyquist criterion
    # Note: B given in Hz. The oscillation frequency
    # is simply the energy.


    E_max = B*Jmax*(Jmax+1);
    if (not for_cos2d):
        E_second = B*(Jmax-2)*(Jmax-1);
    else:
        E_second = B*(Jmax-4)*(Jmax-3);

    frequency = E_max-E_second;
    sampling_frequency = frequency*samples_per_period;
    dt = 1.0/(sampling_frequency);

    return dt;

def xc_window(dasfile):
    pulse_data = np.genfromtxt(dasfile, delimiter=",")
    t = pulse_data[:, 0]
  #  A = np.abs([t[0], t[-1]])

#    if A[0] != A[1]: #if the time axis is not symmetric around 0
#        Q=np.argmax(A)
#        if Q==0: #is the first one the biggest? (so t is negative)
##            t_min=-np.abs(t[0])-1
 #           t_max=np.abs(t[0])+1
 #       else: #or is t mostly positive?
 #           t_min=-np.abs(t[-1])-1
 #           t_max=np.abs(t[-1])+1
#
 #       a=np.array([t_min, 0])
  #      b=np.array([t_max, 0])

   #     pulse_data=np.insert(pulse_data,0,a,axis=0)
    #    pulse_data=np.insert(pulse_data,len(pulse_data),b,axis=0)

      #  np.savetxt(dasfile, pulse_data, delimiter=",") #so we force the file to be symmetric, and edit the input file!!!

    half_range=(np.abs(pulse_data[0,0])+np.abs(pulse_data[-1,0]))/2
    return half_range*1e-12, t[0]*1e-12, t[-1]*1e-12 #return some range, along with the start and end. Maybe change this to just start and end
