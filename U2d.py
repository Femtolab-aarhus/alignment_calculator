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

import interaction
import argparse
import numpy

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate cos^2 theta 2 matrix elements.');

    parser.add_argument('Jmax',type=int,help='J of highest basis vector.');
    parser.add_argument('K',type=int,help='The K quantum number');
    parser.add_argument('M',type=int,help='The M quantum number');

    args = parser.parse_args();

    Jmax = args.Jmax;
    K = int(args.K);
    M = int(args.M);

    if (Jmax < 0):
        raise RuntimeError('J must be positive');

   
    KMsign = numpy.sign(K*M);
    if (KMsign == 0):
        KMsign = 1;

    Kp = abs(K);
    Mp = abs(M);

    if (Kp != 0):
        raise NotImplementedError("Calculation only implemented for linear molecules.");

    interaction.PrepareMeanCos2dMatrix(Jmax,Kp,Mp,KMsign);



