# Personnalize the variable BIN
# Use make to compile the program
#     make clean to remove object files
#     make distclean to remove object files in the current and in the library directory

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
F90FLAGS = -O3 -ip -static
F77 = $(F90)
F77FLAGS = $(F90FLAGS)

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
.PRECIOUS: libAssociate.a
.SUFFIXES: .f .f90 .o .a

.f.o:
	$(F77) $(F77FLAGS) -c $<

.f.a:
	$(F77) $(F77FLAGS) -c $<
	ar Urv $@ $%
	@rm -f $%

.f90.o:
	$(F90) $(F90FLAGS) -I $(MODULES) -I $(MODSTRUCT) -module $(MODULES) -c $<

.f90.a:
	$(F90) $(F90FLAGS) -I $(MODULES) -I $(MODSTRUCT) -module $(MODULES) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	$(BIN) associateConf
associateConf: $(BIN)/associateConf

$(BIN):
	@mkdir $(BIN)

clean:
	@rm -f *.o *.il *.mod $(MODULES)/*.mod libAssociate.a

distclean: clean
	@(cd $(STRUCTDIR)  && $(MAKE) clean)

$(BIN)/associateConf:	$(LIBSTRUCT) libAssociate.a mainAssociateConf.o
	$(LD) $(LDFLAGS) -o $(BIN)/associateConf mainAssociateConf.o libAssociate.a $(LIBSTRUCT)


libAssociate.a: libAssociate.a(data.o) \
		libAssociate.a(lecture.o) \
		libAssociate.a(associate.o)

$(LIBSTRUCT): force_look
	(cd $(STRUCTDIR) && $(MAKE))
