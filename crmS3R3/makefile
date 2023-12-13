# Which targets aren't actual filenames
.PHONY: default test cleanup clean mrproper testaotus

# Set path of implementation of integrator
ifndef INTEGRATOR
#INTEGRATOR = BLieDF
#INTEGRATORP = ../BLieDF
#INTEGRATOR = gena
#INTEGRATORP = ../gena
INTEGRATOR = RATTLie
INTEGRATORP = ../RATTLie
#INTEGRATOR = varint4lie
#INTEGRATORP = ../varint4lie
#INTEGRATOR = SHAKELie
#INTEGRATORP = ../SHAKELie
endif

# Set path of the implementation of the aotus library
AOTP = ../aotus

# Set path of the implementation of the Lie group functions
LIEFUNP = ../liegroup
LIEFUNS = cross_functions.o quaternion_functions.o s3_functions.o s3sdr3_functions.o singular_functions.o

# Set path of the implementation of expandconf
EXPANDCONFIGP = ../expandconfig

# Set Matlab command
MATLAB = matlab -nodesktop

# Set path of the implementation of readLua
READLUAP = ../readLua-for-Matlab-and-Octave

# Fortran compiler:
FC = gfortran
# Flags for Fortran compiler
FFLAGS = -cpp -ffree-line-length-none
FFLAGS := $(FFLAGS) -DINTEGRATOR=$(INTEGRATOR) -DINT_$(INTEGRATOR)
FFLAGS := $(FFLAGS) -DNO_Z
#-DHIGHER_ORDER
#-DLIN_LOG -DAPPROX_TAN_INV_T

# Extra flags for the Fortran compiler (are used for the integrator, as well)
EXTRAFFLAGS := -O3
#EXTRAFFLAGS := -Dpure='' -g -fbounds-check -fimplicit-none -fbacktrace -fcheck=all -finit-real=NaN
ifdef RELAXATION
EXTRAFFLAGS := $(EXTRAFFLAGS) -DRELAXATION
endif
ifdef TWO_ITERATIONS_MINIMUM
EXTRAFFLAGS := $(EXTRAFFLAGS) -DTWO_ITERATIONS_MINIMUM
endif
ifdef VARIABLE_STEPS
EXTRAFFLAGS := $(EXTRAFFLAGS) -DVARIABLE_STEPS
endif
ifdef USE_INDEX_1_APPROX
EXTRAFFLAGS := $(EXTRAFFLAGS) -DUSE_INDEX_1_APPROX
endif
ifdef ONLY_LM_LP
EXTRAFFLAGS := $(EXTRAFFLAGS) -DONLY_LM_LP
endif

#-finit-integer=-77 -DDEBUG_PRINT_ITERATION_MATRIX_AT=5.0
#-DITERATIVE_RELAXATION
#-ffpe-trap=zero,invalid,overflow,underflow
export EXTRAFFLAGS

# Update FFLAGS to use EXTRAFFLAGS
FFLAGS := $(FFLAGS) $(EXTRAFFLAGS)

# Deprecated / TODO remove
#FFLAGS := $(FFLAGS) -Dalt
#FFLAGS := $(FFLAGS) -pg
FFLAGS := $(FFLAGS) -Dzeta1

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

# Linker: (gfortran automatically links the libs required for Fortran)
LD = gfortran
# Flags for Linker
LFLAGS = -llapack
#LFLAGS := $(LFLAGS) -pg

############################################################# TARGETS ##
# default: This target will be built if make is run without arguments
default: crm makefile

# Target to run tests
test: crm test/expandconfig test/als/readLua.m
	cd test && ./start_test.sh

# Target to try if the program works
try: crm test/expandconfig
	cd test && ./start_try.sh

# TODO: Preliminary target
testaotus: obj/testaotus.o obj/libaotus.a
	$(LD) $(LFLAGS) $^ -o test/testaotus
	cd test && ./testaotus

# Target that builds libaotus.a
obj/libaotus.a:
	cd $(AOTP) && ./waf configure build
	cp $(AOTP)/build/libaotus.a  obj/
	cp $(AOTP)/build/*.mod  obj/

# Target that builds $(INTEGRATOR).o and $(INTEGRATOR).mod
obj/$(INTEGRATOR).o:
	$(MAKE) -C $(INTEGRATORP)
	cp $(INTEGRATORP)/obj/$(INTEGRATOR).o   obj/
	cp $(INTEGRATORP)/obj/$(shell echo $(INTEGRATOR) | tr '[:upper:]' '[:lower:]').mod obj/

# Target that builds the Lie group function objects
$(addprefix obj/,$(LIEFUNS)):
	$(MAKE) -C $(LIEFUNP)
	cp $(LIEFUNP)/$@ --target-directory=obj/
	cp $(LIEFUNP)/$(patsubst %.o,%.mod,$@) --target-directory=obj/

# Target that builds get_line_of_variable_length.o
obj/get_line_of_variable_length.o: src/get_line_of_variable_length.F90
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<

# Target that builds main.o
obj/main.o: src/main.F90 obj/libaotus.a obj/get_line_of_variable_length.o
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<

# Target that builds crm.o
obj/crm.o: src/crm.F90 obj/$(INTEGRATOR).o $(addprefix obj/,$(LIEFUNS))
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<

# Target that links everything and builds the executable
crm: obj/$(INTEGRATOR).o obj/crm.o obj/main.o  obj/get_line_of_variable_length.o obj/libaotus.a $(addprefix obj/,$(LIEFUNS))
	$(LD) $(INCLUDE) obj -o $@ $+ $(LFLAGS)

# Target that builds and gets readLua
test/als/readLua.m:
	cd $(READLUAP) && $(MATLAB) -nodesktop -nojvm -r "make,quit"
	cp $(READLUAP)/readLua.m* test/als/

# Target that builds and gets expandconfig
test/expandconfig:
	$(MAKE) -C $(EXPANDCONFIGP)
	cp $(EXPANDCONFIGP)/expandconfig $@

# Target to clean up the object folder
cleanup:
	-rm obj/*

# Target that cleans the executable
clean: cleanup cleanintegrator
	-rm crm
	-rm test/expandconfig

# Target that cleans the $(INTEGRATOR) project folder
cleanintegrator:
	$(MAKE) -C $(INTEGRATORP) clean
	-rm obj/$(INTEGRATOR).o

# Target that cleans only the output folder
cleantest:
	-rm test/out/*
	-rm test/cfg_exp/*

# Target that even almost cleans everything
mrproper: clean cleantest cleanintegrator
	$(MAKE) -C  $(LIEFUNP) clean

# Target that cleans readlua.m
cleanreadlua:
	rm test/als/readLua.m* 2> /dev/null; true
