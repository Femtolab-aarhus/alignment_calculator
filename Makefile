
all:	libpropagation.so libU2d.so  # For Linux/BSD
#all:	libpropagation.dylib libU2d.dylib  # For Darwin (Mac)
#all:	libpropagation.dll libU2d.dll  # For windows

#CC = clang
CC = /usr/local/gfortran/bin/gcc
#CC = i686-w64-mingw32-gcc # for windows, 32 bit
#CC = x86_64-w64-mingw32-gcc # 64 bit

# Fortran compiler to use
FORT = /usr/local/gfortran/bin/gfortran
#FORT = x86_64-w64-mingw32-gfortran # Windows targets, 32/64 bit
#FORT = i686-w64-mingw32-gfortran

WARNINGS = -std=c11 -pedantic -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -D_POSIX_C_SOURCE=20160116

# -fstack-protector-all  is useful for debugging.

DEBUG = #-g

SSE = #-msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -msse4 -msse4a
AVX = #-mavx -mavx2 -mavx512f -mavx512pf -mavx512er -mavx512cd

CFLAGS = -I /Users/adam/anaconda/include/  $(DEBUG) $(WARNINGS) -Wall -mtune=native -march=native $(SSE) $(AVX) -pipe -fpic -Ofast -ffast-math -fassociative-math -funroll-loops -fuse-linker-plugin -frename-registers -fweb -fomit-frame-pointer -funswitch-loops -funsafe-math-optimizations -fno-common

OPENMP = -fopenmp # Comment out to remove openmp support
DISABLE_GSL = # -DNO_GSL # Uncomment to remove code that depends on the GSL.
# Also remove -lgslcblas and -lgsl from LDFLAGS.

BLAS = -lgslcblas
#BLAS = -lblas  # Use another, possibly faster BLAS implementation

LDFLAGS = $(DEBUG) -lgsl $(BLAS) -lm -L /Users/adam/anaconda/lib/



libpropagation.dll:	propagation.o
	@echo 
	@echo Making libpropagation.dll \(For Windows\)
	@echo 
	$(CC) propagation.o $(LDFLAGS) -shared -Wl,--subsystem,console,--out-implib,libpropagation.dll,--add-stdcall-alias -o libpropagation.dll

libpropagation.dylib:	propagation.o
	@echo 
	@echo Making libpropagation.dylib \(For Mac\)
	@echo 
	$(CC) propagation.o $(LDFLAGS) -shared -Wl,-install_name,libpropagation.dylib -static-libgcc -o libpropagation.dylib


libpropagation.so:	propagation.o
	@echo 
	@echo Making libpropagation.so
	@echo 
	$(CC) propagation.o $(LDFLAGS) -shared -Wl,-soname,libpropagation.so -o libpropagation.so

libU2d.dll:	U2d.o wigner/libwigner.a
	@echo 
	@echo Making U2d \(cos^2 theta 2d matrix\) library
	@echo 
	$(CC) -o libU2d.dll -shared U2d.o $(LDFLAGS) -Wl,--subsystem,console,-no-seh,--out-implib,libU2d.dll,--add-stdcall-alias -static-libgcc -L./wigner -lwigner -Wl,--Bstatic -lgfortran -lquadmath $(OPENMP)


libU2d.dylib:	U2d.o wigner/libwigner.a
	@echo 
	@echo Making U2d \(cos^2 theta 2d matrix\) library
	@echo 
	$(CC) U2d.o $(LDFLAGS) -shared -Wl,-install_name,libU2d.dylib  -lwigner -L./wigner -lgfortran -lquadmath -static-libgcc $(OPENMP) -o libU2d.dylib


libU2d.so:	U2d.o wigner/libwigner.a
	@echo 
	@echo Making U2d \(cos^2 theta 2d matrix\) library
	@echo 
	$(CC) U2d.o $(OPENMP) $(LDFLAGS) -shared -Wl,-soname,libU2d.so.1 -lwigner -L./wigner -lgfortran -o libU2d.so




test_propagation: test_propagation.c propagation.c
	cat propagation.c | sed 's/static//g' | sed 's/inline//g' > propagation_nostatic.c
	$(CC) -pg $(DISABLE_GSL) $(filter-out -fomit-frame-pointer, $(CFLAGS)) $(LDFLAGS) -o test_propagation test_propagation.c


propagation.o: propagation.c
	$(CC) $(DISABLE_GSL) $(CFLAGS) -o propagation.o -c propagation.c

propagation.s: propagation.c
	$(CC) $(DISABLE_GSL) $(CFLAGS) -S -fverbose-asm -o propagation.s -c propagation.c

U2d.o: U2d.c
	$(CC) -c $(CFLAGS) $(DISABLE_GSL) $(OPENMP) -o U2d.o U2d.c


clean:
	rm -f *.o libpropagation.so libU2d.so test_propagation propagation_nostatic.c
	(cd wigner && $(MAKE) clean)


wigner/libwigner.a: $(wildcard wigner/*.f)
	(cd wigner && $(MAKE) CC=$(CC) COMPILER=$(FORT))


