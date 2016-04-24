#CC = clang
CC = gcc

#WINCC = i686-w64-mingw32-gcc
WINCC = x86_64-w64-mingw32-gcc

WARNINGS = -std=c11 -pedantic -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -D_POSIX_C_SOURCE=20160116

CFLAGS = $(WARNINGS) -Wall -mtune=native -march=native -pipe -fpic -Ofast -ffast-math -fassociative-math -funroll-loops -fuse-linker-plugin -frename-registers -fweb -fomit-frame-pointer -funswitch-loops -funsafe-math-optimizations -fno-common
LDFLAGS = -lm 


all:	libpropagation.so libU2d.so

libpropagation.dll:	propagation_win.o
	@echo 
	@echo Making libpropagation.dll \(windows\)
	@echo 
	$(WINCC) $(LDFLAGS) -shared -Wl,--subsystem,windows,--out-implib,libpropagation.dll,--add-stdcall-alias propagation_win.o -o libpropagation.dll


libpropagation.so:	propagation.o
	@echo 
	@echo Making libpropagation.so
	@echo 
	$(CC) $(LDFLAGS) -shared -Wl,-soname,libpropagation.so.1 propagation.o -o libpropagation.so



propagation.o: propagation.c
	$(CC) $(CFLAGS) -o propagation.o -c propagation.c

propagation.s: propagation.c
	$(CC) $(CFLAGS) -S -fverbose-asm -o propagation.s -c propagation.c

propagation_win.o: propagation.c
	$(WINCC) $(CFLAGS) -o propagation_win.o -c propagation.c

U2d.o: U2d.c
	$(CC) -c $(CFLAGS) -fopenmp -o U2d.o U2d.c
U2d_win.o: U2d.c
	$(WINCC) -c $(CFLAGS) -DNO_GSL -fopenmp -o U2d_win.o U2d.c


clean:
	rm -f *.o libpropagation.so libU2d.so
	(cd wigner && make clean)


wigner/libwigner.a: $(wildcard wigner/*.f)
	(cd wigner && make)

libU2d.so:	U2d.o wigner/libwigner.a
	@echo 
	@echo Making U2d \(cos^2 theta 2d matrix\) library
	@echo 
	$(CC) $(LDFLAGS) -shared -Wl,-soname,libU2d.so.1 U2d.o -lgfortranbegin -lgfortran -lwigner -L./wigner -lgsl -lgslcblas -fopenmp -o libU2d.so

libU2d.dll:	U2d_win.o wigner/libwigner.a
	@echo 
	@echo Making U2d \(cos^2 theta 2d matrix\) library
	@echo 
	$(WINCC) $(LDFLAGS) -shared -Wl,--subsystem,windows,--out-implib,libU2d.dll,--add-stdcall-alias U2d_win.o -lwigner -lgfortran -lquadmath -L./wigner -fopenmp -o libU2d.dll



