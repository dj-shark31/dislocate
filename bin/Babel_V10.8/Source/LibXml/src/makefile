.PRECIOUS: libXmlParse.a
.SUFFIXES: .f90 .o .a

MODULES = ./Modules

## #F90	=	f90
## #F90	=	lf95
## F90	=	ifort
## #F90FLAGS	=	-c
## #F90	=	gfortran
## F90FLAGS	=	-c


.f90.o:
	$(F90) $(F90FLAGS) -I $(MODULES) -c $<

.f90.a:
	$(F90) $(F90FLAGS) -I $(MODULES) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	libXmlParse.a 

clean:	
	@rm -f *.o *.mod $(MODULES)/*.mod libXmlParse.a

libXmlParse.a:	libXmlParse.a(xmlparse.o) \
       libXmlParse.a(read_xml_prims.o)

tstparse:	tstparse.o libXmlParse.a
	$(F90) -o tstparse tstparse.o libXmlParse.a

xmlreader:	xmlreader.o libXmlParse.a
	$(F90) -o xmlreader xmlreader.o libXmlParse.a

