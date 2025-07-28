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

##################################
.SUFFIXES: .f .f90 .o .a

.f.o:
	$(F77) $(F77FLAGS) -c $<

.f.a:
	$(F77) $(F77FLAGS) -c $<
	ar Urv $@ $%
	@rm -f $%

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f90.a:
	$(F90) $(F90FLAGS) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	vitek_bcc vitek_bcc_extract vitek_bcc_M111 vitek_hcp vitek_hcp_extract
vitek_bcc:	$(BIN)/vitek_bcc
vitek_bcc_extract:	$(BIN)/vitek_bcc_extract
vitek_bcc_M111:	$(BIN)/vitek_bcc_M111
vitek_hcp:	$(BIN)/vitek_hcp
vitek_hcp_extract:	$(BIN)/vitek_hcp_extract

$(BIN):
	@mkdir $(BIN)

clean:
	@rm -f *.o 

$(BIN)/vitek_bcc:	mainVitek_bcc.o
	$(LD) $(LDFLAGS) -o $(BIN)/vitek_bcc mainVitek_bcc.o

$(BIN)/vitek_bcc_extract:	mainVitek_bcc_extract.o
	$(LD) $(LDFLAGS) -o $(BIN)/vitek_bcc_extract mainVitek_bcc_extract.o

$(BIN)/vitek_bcc_M111:	mainVitek_bcc_M111.o
	$(LD) $(LDFLAGS) -o $(BIN)/vitek_bcc_M111 mainVitek_bcc_M111.o

$(BIN)/vitek_hcp:	mainVitek_hcp.o
	$(LD) $(LDFLAGS) -o $(BIN)/vitek_hcp mainVitek_hcp.o

$(BIN)/vitek_hcp_extract:	mainVitek_hcp_extract.o
	$(LD) $(LDFLAGS) -o $(BIN)/vitek_hcp_extract mainVitek_hcp_extract.o

