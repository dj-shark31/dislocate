# Personnalize the variable BIN
# Use make to compile the program
#     make clean to remove object files
#     make distclean /make veryclean to remove object files in the current and in the library directory

BIN = ../../../bin
MODULES = ./Modules

# Intel fortran compiler
F90 = ifort
# Debugging on Linux (ifort)
#F90FLAGS = -O0 -fpe0 -CB -check all -e95 -std95 -traceback -warn all
#F90FLAGS = -O0 -fpe0 -CB -e95 -std95 -traceback -warn all
#F90FLAGS = -O0 -fpe0 -CB -traceback -warn all
# Optimizing on Linux (ifort)
#  Do not use option -axN
F90FLAGS = -O3 -ip 
F77 = $(F90)
F77FLAGS = $(F90FLAGS)

### End of options #########################################################

STRUCTDIR = ../../LibStruct
LIBSTRUCT = $(STRUCTDIR)/libStruct.a
MODSTRUCT = $(STRUCTDIR)/Modules

##################################
LD = $(F90)
LDFLAGS = $(F90FLAGS)
export F90
export F90FLAGS
export F77
export F77FLAGS

##################################
.PHONY: force_look
.PRECIOUS: libBabel.a
.SUFFIXES: .f .f90 .o .a

.f.o:
	$(F77) $(F77FLAGS) -module $(MODULES) -c $<

.f.a:
	$(F77) $(F77FLAGS) -module $(MODULES) -c $<
	ar Urv $@ $%
	@rm -f $%

.f90.o:
	$(F90) $(F90FLAGS) -I $(MODULES) -I $(MODSTRUCT) -module $(MODULES) -c $<

.f90.a:
	$(F90) $(F90FLAGS) -I $(MODULES) -I $(MODSTRUCT) -module $(MODULES) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	$(BIN) buildXyzLine
buildXyzLine:	$(BIN)/buildXyzLine

$(BIN):
	@mkdir $(BIN)
clean:
	@rm -f *.o *.il *.mod $(MODULES)/*.mod lib*.a

veryclean: distclean
distclean: clean
	@(cd $(STRUCTDIR)  && $(MAKE) clean)

$(BIN)/buildXyzLine:	$(LIBSTRUCT) mainBuildXyzLine.o
	$(LD) $(LDFLAGS) -o $(BIN)/buildXyzLine mainBuildXyzLine.o $(LIBSTRUCT)

$(LIBSTRUCT): force_look
	(cd $(STRUCTDIR) && $(MAKE))
