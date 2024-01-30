
EXE_DIR := bin/

PROG:=$(EXE_DIR)a.out

FC := gfortran
FFLAGS := -fconvert=big-endian -ffree-line-length-none #-mcmodel=large 
FFLAGS += -Wall -fbounds-check -O  -fbacktrace -g  -Wuninitialized #-ffpe-trap=invalid,zero,overflow,underflow

#FC:=/data/home/sfujibayashi/libs/hdf5-intel/bin/h5fc
# FC:=ifort
# FFLAGS0:= -convert big_endian -mcmodel=large -shared-intel -fpic # -qopenmp
# FFLAGS := $(FFLAGS0) -O3 -xHost  -O0 -CB -traceback -g -fpe0 

.SUFFIXES: .f90 .F90 .o

SRC :=  module_ptf_reaclib.f90 polint.f90 locate.f90 const_mod.f90 module_nse.f90 main.f90

SRC_DIR := src/

OBJ_DIR := obj/

FFLAGS += -J$(OBJ_DIR)

SRC := $(addprefix $(SRC_DIR), $(SRC))

OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC:.f90=.o)))
# MOD_FILES := $(addprefix $(OBJ_DIR),$(notdir $(MOD:.f90=.mod)))




.PHONY: all test

all: dirs $(PROG)

test:
	@echo $(SRC)
	@echo $(OBJ_FILES)


dirs : $(EXE_DIR) $(OBJ_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(PROG): $(OBJ_FILES)
	$(FC) $(FFLAGS) -o $@  $(OBJ_FILES) $(OMPFLAGS) $(LIBS)

$(OBJ_DIR)%.o: $(SRC_DIR)%.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)%.mod: $(SRC_DIR)%.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) $(LIBS) -c $< -o $@


# clean rule
.PHONY: clean
clean:
	$(RM) $(OBJ_FILES) fort* *.mod
