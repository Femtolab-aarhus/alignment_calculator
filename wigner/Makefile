COMPILER = gfortran
#COMPILER = x86_64-w64-mingw32-gfortran
#COMPILER = i686-w64-mingw32-gfortran
#COMPILER = ifort
CC = gcc
#CC = icc
#CC = x86_64-w64-mingw32-gcc
#CC = i686-w64-mingw32-gcc

CHECK = -fcheck=array-temps,bounds #-ftree-vectorizer-verbose=1
# 
FFLAGS = -fPIC -O3 -fopenmp -fexternal-blas -Wall -pedantic -mtune=native -march=native -pipe -ffast-math -fuse-linker-plugin -fweb -frename-registers -funroll-loops -funroll-all-loops $(CHECK)
# -g and -flto are not compatible
#FFLAGS = -g -O0 -fopenmp -Wall -fexternal-blas -pedantic -mtune=native -march=native -pipe -fcheck=all 

INCLUDE = -I/usr/include
LINK = -L/usr/lib/ -lm -llapack -lmpi -lgsl -lgslcblas -lblas -llapack

CFLAGS = -fPIC -O3 -Wall -pedantic -mtune=native -march=native -pipe -ffast-math -fuse-linker-plugin -fweb -frename-registers -std=gnu99 -funroll-loops
#CFLAGS = -g -O0 -Wall -pedantic -mtune=native -march=native -pipe -ffast-math -std=gnu99

source = $(wildcard *.f)
source90 = $(wildcard *.f90)
obj = $(source:.f=.o) $(source90:.f90=.o)

csource = $(wildcard *.c)
cobj = $(csource:.c=.o)

all: libwigner.a


libwigner.a: $(obj) $(cobj)
	ar rcs libwigner.a $(obj) $(cobj)

%.o : %.f
	$(COMPILER) $(FFLAGS) $(INCLUDE) -o $@ -c $<

%.o : %.f90
	$(COMPILER) $(FFLAGS) $(INCLUDE) -o $@ -c $<

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

.PHONY: clean
clean:
	rm -f *.o libwigner.a
