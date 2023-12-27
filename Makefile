
FC = mpif90 #-fc=gfortran-7.2.0
#FC = gfortran-7.2.0
#MPIDIR = /usr/pack/mpich-3.1.1-ma
#-I$(MPIDIR)/include
FCOPT = -O3 -fopenmp 
LNFLAG = -fopenmp
#LIBS = -Wl,-rpath,$(MPIDIR)/lib  -L$(MPIDIR)/lib -lmpi -lmpifort 

#FCOPT = -g -O0 -fbounds-check -fmax-errors=3
#LIBS = -L../lapack-3.8.0 -llapack -lrefblas
SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)

TARGET = poisson_p

all: $(TARGET)


clean:
	rm *.o *.mod $(TARGET)

mpiclean:
	rm out.* err.*

poisson.o: cg.o matdef.o adj_map.o precision.o sparsealg.o
poisson.o: mpi_globals.o clock.o
matdef.o: precision.o
sparsealg.o: precision.o matdef.o
cg.o: precision.o matdef.o sparsealg.o mpi_globals.o
#mpi_globals.o : mpi.o

$(TARGET): $(OBJS)
	$(FC) $(LNFLAG)  -o $(TARGET) $(OBJS) $(LIBS)


%.o: %.f90
	$(FC) -c $(FCOPT) $<




