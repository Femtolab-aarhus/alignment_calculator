
all:	libpropagation.so libU2d.so  # For POSIX
#all:	libpropagation.dll libU2d.dll  # For windows

#CC = clang
CC = gcc
#CC = i686-w64-mingw32-gcc # for windows, 32 bit
#CC = x86_64-w64-mingw32-gcc # 64 bit

# Fortran compiler to use
FORT = gfortran
#FORT = x86_64-w64-mingw32-gfortran # Windows targets, 32/64 bit
#FORT = i686-w64-mingw32-gfortran

WARNINGS = -std=c11 -pedantic -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -D_POSIX_C_SOURCE=20160116

# -fstack-protector-all  is useful for debugging.

CFLAGS = $(WARNINGS) -Wall -mtune=native -march=native -msse2 -pipe -fpic -Ofast -ffast-math -fassociative-math -funroll-loops -fuse-linker-plugin -frename-registers -fweb -fomit-frame-pointer -funswitch-loops -funsafe-math-optimizations -fno-common

OPENMP = -fopenmp # Comment out to remove openmp support
DISABLE_GSL = # -DNO_GSL # Uncomment to remove code that depends on the GSL.
# Also remove -lgslcblas and -lgsl from LDFLAGS.

LDFLAGS = -lm -lgslcblas -lgsl 



libpropagation.dll:	propagation.o
	@echo 
	@echo Making libpropagation.dll \(For Windows\)
	@echo 
	$(CC) propagation.o $(LDFLAGS) -shared -Wl,--subsystem,console,--out-implib,libpropagation.dll,--add-stdcall-alias -o libpropagation.dll


libpropagation.so:	propagation.o
	@echo 
	@echo Making libpropagation.so
	@echo 
	$(CC) propagation.o $(LDFLAGS) -shared -Wl,-soname,libpropagation.so.1 -o libpropagation.so

libU2d.dll:	U2d.o wigner/libwigner.a
	@echo 
	@echo Making U2d \(cos^2 theta 2d matrix\) library
	@echo 
	$(CC) -o libU2d.dll -shared U2d.o $(LDFLAGS) -Wl,--subsystem,console,-no-seh,--out-implib,libU2d.dll,--add-stdcall-alias -static-libgcc -L./wigner -lwigner -Wl,--Bstatic -lgfortran -lquadmath $(OPENMP)

libU2d.so:	U2d.o wigner/libwigner.a
	@echo 
	@echo Making U2d \(cos^2 theta 2d matrix\) library
	@echo 
	$(CC) $(LDFLAGS) -shared -Wl,-soname,libU2d.so.1 U2d.o -lgfortranbegin -lgfortran -lwigner -L./wigner $(OPENMP) -o libU2d.so





propagation.o: propagation.c
	$(CC) $(DISABLE_GSL) $(CFLAGS) -o propagation.o -c propagation.c

propagation.s: propagation.c
	$(CC) $(DISABLE_GSL) $(CFLAGS) -S -fverbose-asm -o propagation.s -c propagation.c

U2d.o: U2d.c
	$(CC) -c $(CFLAGS) $(DISABLE_GSL) $(OPENMP) -o U2d.o U2d.c


clean:
	rm -f *.o libpropagation.so libU2d.so
	(cd wigner && $(MAKE) clean)


wigner/libwigner.a: $(wildcard wigner/*.f)
	(cd wigner && $(MAKE) CC=$(CC) COMPILER=$(FORT))


