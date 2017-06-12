#copied from http://www.<webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
# A simple hand-made makefile for a package including applications
# built from Fortran 90 sources, taking into account the usual
# dependency cases.

# This makefile works with the GNU make command, the one find on
# GNU/Linux systems and often called gmake on non-GNU systems, if you
# are using an old style make command, please see the file
# Makefile_oldstyle provided with the package.

# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check -ffree-form -fbacktrace -Wall -Wtabs -fcheck=pointer
#FCFLAGS += -O2
#flags for avoiding errors
#FCFLAGS += -pedantic-errors

# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I.
# libraries needed for linking, unused in the examples
# modules here /ingrid
LDFLAGS = -llapack -lblas

#precompiled modules here. E.g. common.
#MODULES =$(wildcard ../common/*.o)
#MODULES = module1.o

# List of executables to be built within the package
PROGRAMS = main

#list of objects to be compiled
# "make" builds all
OBJECTS = module1.o main.o 


all: $(PROGRAMS)


$(PROGRAMS):$(OBJECTS)
	$(FC) $(FCFLAGS) $^ -o $@  $(MODULES)  $(LDFLAGS)


# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

# Utility targets
#.PHONY: clean veryclean
.PHONY:clean

clean:
	rm -f *.o *.mod *.MOD

#veryclean: clean
#	rm -f *~ $(PROGRAMS)
