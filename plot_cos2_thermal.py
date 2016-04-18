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

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import config
import numpy 
from numpy import pi,exp,sqrt,log,sin,cos
import argparse
import time
import datetime

parser = argparse.ArgumentParser(description='Plot cos^2 trace of molecule(s) calculated with cos2_thermal.py.');

parser.add_argument('plot_title',nargs=1,type=str,help='Title of plot.');
parser.add_argument('file',nargs="+",type=str,help='File(s) to plot. Add :legend to plot it with that legend.');
parser.add_argument('--dontshow',action="store_true",help='Just save the plot, don\'t show it.');


args = parser.parse_args();
title = args.plot_title[0];



size = numpy.array([18.5,10.5])*0.8;

fig = plt.figure(figsize=size, dpi=80);
plt.ylim([0,1]);
#plt.ylim([0.25,0.9]);
#plt.xlim([75,79]);

#plt.axes((0.072,0.1,0.9,0.85));

matplotlib.rcParams.update({'font.size': 24})

print(args.file);

colors = ['b','r','g','m','c','k'];

c = 0;
for filename in args.file:
    
    print(filename);
    filename = filename.split(":");
    if len(filename) > 1:
        legend = filename[1];
    else:
        legend = None;
    filename = filename[0];

    print("Loading " + filename);
    contents = numpy.load(filename);

    t = contents["t"];
    cos2 = contents["cos2"];
    meta = contents["meta"];

    plt.plot(t*1e12,cos2,label=legend,color=colors[c%len(colors)]);
    #plt.plot(t*1e12,cos2,'o',color=colors[c%len(colors)]);
    c = c + 1;

plt.xlabel('Time [ps]');
plt.ylabel('<cos^2 theta>');
plt.legend(fontsize=20,frameon=False);
plt.grid();


plt.title(title);

print("Saving image to plot.pdf");
plt.savefig("plot.pdf",bbox_inches='tight');

if (not args.dontshow):
    plt.show();

