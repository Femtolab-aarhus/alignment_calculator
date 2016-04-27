#!/bin/bash

# It may be a good idea to run this script before simulations.
# It disables threading in common libraries, so we avoid running
# more than 1 thread per computation unit.
# For this script to work, you must run it with "source".
# I.e:
# source env.sh

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

