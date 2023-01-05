PROG=a.out

# FC = gfortran
# FFLAGS0= -fconvert=big-endian -ffree-line-length-none #-mcmodel=large 
# FFLAGS = -Wall -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace -g $(FFLAGS0)

#FC=/data/home/sfujibayashi/libs/hdf5-intel/bin/h5fc
FC=ifort
FFLAGS0= -convert big_endian -mcmodel=large -shared-intel -fpic # -qopenmp
FFLAGS = $(FFLAGS0) -O3 -xHost  -O0 -CB -traceback -g -fpe0 

.SUFFIXES: .f90 .F90 .o

#SRC1 = eostab_mod.f90 unit_mod.f90 const_mod.f90 main.f90 read_helm_table.f90 helmeos.f90 timeos.f90
#SRC1 = const_mod.f90 readmintem.f90 timeos.f90
#SRC1 = const_mod.f90 read.f90 timeos.f90
#SRC1 = const_mod.f90 readbeta.f90 timeos.f90
SRC1 = const_mod.f90 module_nse.f90 main.f90
#SRC1 = EOSformat.f90


OBJS = $(SRC1:%.f90=%.o)

.PHONY: all
all: $(PROG)

# 1st rule
$(PROG): $(OBJS)
	$(FC) $(FFLAGS) -o $(PROG) $(OMPFLAGS) $(LIBS) $(OBJS)

.f90.o:
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $<

.c.o:
	$(FCC) -O -c $<

.f90.s:
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -S $<

.F90.o:
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $<

.F90.s:
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -S $<

beta.o: beta.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $<

module_nuclei_quantum.o: module_nuclei_quantum.F90
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $<
module_nse.o: module_nse.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $<
module_sekig_eos.o: module_sekig_eos.F90
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $<


# clean rule
.PHONY: clean
clean:
	$(RM) $(OBJS) fort* *.mod
