.PRECIOUS: libStruct.a
.SUFFIXES: .f .f90 .o .a

MODULES = ./Modules

.f.o:
	$(F77) $(F77FLAGS) -c $<

.f.a:
	$(F77) $(F77FLAGS) -c $<
	ar Urv $@ $%
	@rm -f $%

.f90.o:
	$(F90) $(F90FLAGS) -I $(MODULES) -c $<

.f90.a:
	$(F90) $(F90FLAGS) -I $(MODULES) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	libStruct.a

clean:
	@rm -f *.o *.mod $(MODULES)/*.mod libStruct.a 

libStruct.a: libStruct.a(mathStruct.o) \
	        libStruct.a(comment.o) \
		libStruct.a(xyz.o) \
		libStruct.a(onlyAtoms.o) \
		libStruct.a(cfg.o) \
		libStruct.a(gin.o) \
		libStruct.a(siesta.o) \
		libStruct.a(lisa.o) \
		libStruct.a(ndm.o) \
		libStruct.a(lammps.o) \
		libStruct.a(poscar.o) \
		libStruct.a(transformation.o) \
		libStruct.a(symmetry3D.o) \
		libStruct.a(structure.o)
