# We want to control all multi processing. Use of multiple cores should be managed as far away from tight loops as possible
# We don't want more processes than cores, and we don't want to
# do threading and inter-process communcation inside tight loops.
# That would also hurt CPU caching in several ways.
import os
os.environ["OMP_NUM_THREADS"] = "1";
os.environ["OPENBLAS_NUM_THREADS"]= "1";
os.environ["MKL_NUM_THREADS"] = "1";

# Note that these environment variables typically have to be set before
# the library in question is even loaded. 
