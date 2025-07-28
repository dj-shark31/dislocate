BIN = ../../../bin

# Intel fortran compiler
F90 = ifort
# Debugging on Linux (ifort)
#F90FLAGS = -O0 -fpe0 -CB -check all -e95 -std95 -traceback -warn all
#F90FLAGS = -O0 -fpe0 -CB -e95 -std95 -traceback -warn all
#F90FLAGS = -O0 -fpe0 -CB -traceback -warn all
# Optimizing on Linux (ifort)
#  Do not use option -axN
F90FLAGS = -O3 -ip -static 
F77 = $(F90)
F77FLAGS = $(F90FLAGS)

##################################
LD = $(F90)
LDFLAGS = $(F90FLAGS)
export F90
export F90FLAGS
export F77
export F77FLAGS

STRUCTDIR = ../../LibStruct
LIBSTRUCT = $(STRUCTDIR)/libStruct.a
MODSTRUCT = $(STRUCTDIR)/Modules


##################################
.PHONY: force_look
.SUFFIXES: .f .f90 .o .a

.f.o:
	$(F77) $(F77FLAGS) -c $<

.f.a:
	$(F77) $(F77FLAGS) -c $<
	ar Urv $@ $%
	@rm -f $%

.f90.o:
	$(F90) $(F90FLAGS) -module $(MODSTRUCT) -c $<

.f90.a:
	$(F90) $(F90FLAGS) -module $(MODSTRUCT) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	distances
distances:	$(BIN)/distances

$(BIN):
	@mkdir $(BIN)

clean:
	@rm -f *.o 

$(LIBSTRUCT): force_look
	(cd $(STRUCTDIR) && $(MAKE))

$(BIN)/distances:	$(LIBSTRUCT) distances_data.o mainDistances.o
	$(LD) $(LDFLAGS) -o $(BIN)/distances mainDistances.o distances_data.o $(LIBSTRUCT)

