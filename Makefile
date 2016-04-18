#CC = clang
CC = gcc

#WINCC = i686-w64-mingw32-gcc
WINCC = x86_64-w64-mingw32-gcc

WARNINGS = -std=c11 -pedantic -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -D_POSIX_C_SOURCE=20160116

CFLAGS = $(WARNINGS) -Wall -mtune=native -march=native -pipe -fpic -Ofast -ffast-math -fassociative-math -funroll-loops -fuse-linker-plugin -frename-registers -fweb -fomit-frame-pointer -funswitch-loops -funsafe-math-optimizations -fno-common
LDFLAGS = -lm 


all:	libpropagation.so U2d

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


clean:
	rm -f propagation.o propagation_win.o libpropagation.so


U2d:	U2d.c
	@echo 
	@echo Making U2d program
	@echo 
	$(CC) U2d.c $(CFLAGS) -O3 -lm -lgsl -lgslcblas -fopenmp -o U2d

U2d.exe:	U2d.c
	@echo 
	@echo Making U2d program
	@echo 
	$(WINCC) U2d.c -DWINDOWS $(CFLAGS) -O3 -lm -lgsl -lgslcblas -fopenmp -o U2d.exe



