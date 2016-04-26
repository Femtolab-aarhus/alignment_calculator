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

import argparse
import multiprocessing
# Tell U2dcalc not to limit the number of threads to use
os.environ["OMP_NUM_THREADS"] = str(multiprocessing.cpu_count());
import U2dcalc


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Precalculate cos^2 theta 2d matrix elements for quick lookup instead of online calculation.');

    parser.add_argument('Jmax',type=int,help='Maximum J quantum number.');
    parser.add_argument('Kmax',type=int,help='Maximum K quantum number.');
    parser.add_argument('Mmax',type=int,help='Maximum M quantum number.');


    args = parser.parse_args();

    Jmax = args.Jmax;
    Kmax = args.Kmax;
    Mmax = args.Mmax;

    U2dcalc.precalculate_matrix_elements(Jmax,Kmax,Mmax);


