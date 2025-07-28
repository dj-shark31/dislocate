# Personnalize the variable BIN
# Use make to compile the program
#     make clean to remove object files

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

LD = $(F90)
LDFLAGS = $(F90FLAGS)
export F90
export F90FLAGS
export F77
export F77FLAGS

##################################
.PHONY: force_look
.PRECIOUS: 
.SUFFIXES: .f .f90 .o .a

.f.o:
	$(F77) $(F77FLAGS) -module $(MODULES) -c $<

.f.a:
	$(F77) $(F77FLAGS) -module $(MODULES) -c $<
	ar Urv $@ $%
	@rm -f $%

.f90.o:
	$(F90) $(F90FLAGS) -I $(MODULES) -module $(MODULES) -c $<

.f90.a:
	$(F90) $(F90FLAGS) -I $(MODULES) -module $(MODULES) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	$(BIN) burgers1D burgers2D 
burgers1D:	$(BIN)/burgers1D
burgers2D:	$(BIN)/burgers2D

$(BIN):
	@mkdir $(BIN)

clean:
	@rm -f *.o *.il *.mod $(MODULES)/*.mod lib*.a

veryclean: distclean
distclean: clean

$(BIN)/burgers1D:	mainBurgers1D.o
	$(LD) $(LDFLAGS) -o $(BIN)/burgers1D mainBurgers1D.o

$(BIN)/burgers2D:	mainBurgers2D.o
	$(LD) $(LDFLAGS) -o $(BIN)/burgers2D mainBurgers2D.o
