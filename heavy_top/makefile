# Which targets aren't actual filenames
.PHONY: default test disstest cleanup clean mrproper testaotus try gena BLieDF RATTLie

# Set path of implementation of integrator
# if one of the goals gena, BLieDF or RATTLie is given, the variables for the
# integrator and its path are set automatically
ifneq ($(filter gena,$(MAKECMDGOALS)),)
  INTEGRATOR = gena
  INTEGRATORP = ../gena
else ifneq ($(filter BLieDF,$(MAKECMDGOALS)),)
  INTEGRATOR = BLieDF
  INTEGRATORP = ../BLieDF
else ifneq ($(filter RATTLie,$(MAKECMDGOALS)),)
  INTEGRATOR = RATTLie
  INTEGRATORP = ../RATTLie
else
  # if none of those special goals is given, the integrator can be set here
  # manually
  #INTEGRATOR = BLieDF
  #INTEGRATORP = ../BLieDF
  #INTEGRATOR = gena
  #INTEGRATORP = ../gena
  INTEGRATOR = RATTLie
  INTEGRATORP = ../RATTLie
endif

# Set path of the implementation of the aotus library
AOTP = ../aotus

# Set path of the implementation of the Lie group functions
LIEFUNP = ../liegroup
LIEFUNS = cross_functions.o quaternion_functions.o s3_functions.o s3sdr3_functions.o singular_functions.o

# Set path of the implementation of expandconf
EXPANDCONFIGP = ../expandconfig

# Set Matlab command
MATLAB = matlabnoterm

# Set path of the implementation of readLua
READLUAP = ../readLua-for-Matlab-and-Octave

# Fortran compiler:
FC = gfortran
# Flags for Fortran compiler
FFLAGS = -cpp -ffree-line-length-none
FFLAGS := $(FFLAGS) -DNO_Z -DINTEGRATOR=$(INTEGRATOR) -DINT_$(INTEGRATOR)

# Extra flags for the Fortran compiler (are used for integrator, as well)
EXTRAFFLAGS := -O3
#EXTRAFFLAGS := -Wall -Wextra -Dpure='' -g -fbounds-check -fimplicit-none -fbacktrace -fcheck=all -finit-real=NaN -finit-integer=-77 -ffpe-trap=zero,invalid
#EXTRAFFLAGS := $(EXTRAFFLAGS) -pg
export EXTRAFFLAGS

# Update FFLAGS to use EXTRAFFLAGS
FFLAGS := $(FFLAGS) $(EXTRAFFLAGS)

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
default: heavy_top makefile

# Target to run tests
test: heavy_top test/expandconfig test/als/readLua.m
	cd test && ./start_test.sh

# Target to run tests for the dissertation
disstest: heavy_top test/expandconfig test/als/readLua.m
	cd test && ./start_diss_test.sh

# Target to try if the program works
try:  heavy_top test/expandconfig
	cd test && ./start_try.sh

# Target that builds the Lie group function objects
$(addprefix obj/,$(LIEFUNS)):
	$(MAKE) -C $(LIEFUNP)
	cp $(LIEFUNP)/$@ --target-directory=obj/
	cp $(LIEFUNP)/$(patsubst %.o,%.mod,$@) --target-directory=obj/

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

# Target that builds get_line_of_variable_length.o
obj/get_line_of_variable_length.o: src/get_line_of_variable_length.F90
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<

# Target that builds main.o
obj/main.o: src/main.F90 obj/libaotus.a obj/get_line_of_variable_length.o
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<

# Target that builds lie_group_functions.o
obj/lie_group_functions.o: src/lie_group_functions.F90 $(addprefix obj/,$(LIEFUNS))
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<

# Target that builds heavy_top.o
obj/heavy_top.o: src/heavy_top.F90 obj/$(INTEGRATOR).o obj/lie_group_functions.o $(addprefix obj/,$(LIEFUNS))
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<

# Target that links everything and builds the executable
heavy_top: obj/$(INTEGRATOR).o obj/heavy_top.o obj/lie_group_functions.o obj/main.o obj/get_line_of_variable_length.o obj/libaotus.a $(addprefix obj/,$(LIEFUNS))
	$(LD) $(INCLUDE) obj -o $@ $+ $(LFLAGS)

# Target that builds and gets readLua
test/als/readLua.m:
	#cd $(READLUAP) && $(MATLAB) -nodesktop -nojvm -r "make,quit"
	cp $(READLUAP)/readLua.m* test/als/

# Target that builds and gets expandconfig
test/expandconfig:
	$(MAKE) -C $(EXPANDCONFIGP)
	cp $(EXPANDCONFIGP)/expandconfig $@

# Target to clean up the object folder
cleanup:
	-rm obj/*

# Target that cleans the $(INTEGRATOR) project folder
cleanintegrator:
	$(MAKE) -C $(INTEGRATORP) clean

# Target that cleans the executable
clean: cleanup cleanintegrator
	-rm heavy_top
	-rm test/expandconfig

# Target that cleans only the output folder
cleantest:
	-rm test/out/*
	-rm test/cfg_exp/*

# Target that even cleans everything
mrproper: clean cleantest
	-rm test/als/readLua.m*
