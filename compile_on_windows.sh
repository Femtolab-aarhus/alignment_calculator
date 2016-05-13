#!/bin/sh


#Due to licensing issues with my compiler, I am not compiling in openmp.
#You can off course compile it in on your own.

make libpropagation.dll libU2d.dll CC=x86_64-w64-mingw32-gcc FORT=x86_64-w64-mingw32-gfortran OPENMP=""

