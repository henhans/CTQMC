all: ctqmc

#LIBS = -L$LIBRARY_PATH -I$INCLUDE -lmkl_blas95 -lmkl_intel -lmkl_lapack95 -lmkl_intel_thread -lmkl_core -lm  -liomp5 -lgsl
LIBS = -llapack -lblas -lgsl
CC = icpc#g++
CFLAGS = #-I/usr/include/i386-linux-gnu/c++/4.8/ -Os -g -funroll-loops -Wall -pedantic -D_OMP -openmp

ctqmc: ctqmc.o
	$(CC) $(CFLAGS) ctqmc.o -o ctqmc $(LIBS)

ctqmc.o: ctqmc.cc
	$(CC) $(CFLAGS) -c ctqmc.cc

clean: 
	rm -rf *.o ctqmc

