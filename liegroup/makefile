.PHONY: default clean .FORCE


OBJS = cross_functions.o quaternion_functions.o s3_functions.o s3sdr3_functions.o singular_functions_1e-6.o

# Set Fortran compiler
FC = gfortran

# Set flags used by Fortran
ifeq "$(FC)" "gfortran"
	INCLUDE = -I
	MODULEP = -J
else
ifeq "$(FC)" "ifort"
	INCLUDE = -I
	MODULEP = -module
endif
endif

# Set Fortran flags
FFLAGS = -O3 -Wall -cpp -ffree-line-length-none

# Use exported Fortran flags
ifdef EXTRAFFLAGS
FFLAGS := $(FFLAGS) $(EXTRAFFLAGS)
endif

default: $(addprefix obj/,$(OBJS)) makefile

test:	obj/test.o $(addprefix obj/,$(OBJS)) .FORCE
	$(FC) $(MODULEP) obj $(FFLAGS) -o $@ obj/test.o $(addprefix obj/,$(OBJS))
	./test

clean:
	-rm obj/*.o obj/*.mod

obj/%.o: src/%.F90
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<

# Dependecies
obj/test.o: $(addprefix obj/,$(OBJS))

obj/s3_functions.o: obj/cross_functions.o obj/quaternion_functions.o obj/singular_functions.o

obj/s3sdr3_functions.o: obj/s3_functions.o
