.PRECIOUS: libLapack.a
.SUFFIXES: .f .o .a
MODULES = ./Modules


.f.o:
	$(F77) $(F77FLAGS) -c $<

.f.a:
	$(F77) $(F77FLAGS) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	libLapack.a

clean:
	@rm -f *.o *.mod $(MODULES)/*.mod libLapack.a

libLapack.a: libLapack.a(dgesv.o)\
	libLapack.a(dgetf2.o)\
	libLapack.a(dgetrf.o)\
	libLapack.a(dgetrs.o)\
	libLapack.a(dlabad.o)\
	libLapack.a(dladiv.o)\
	libLapack.a(dlamck.o)\
	libLapack.a(dlapy2.o)\
	libLapack.a(dlapy3.o)\
	libLapack.a(dlaswp.o)\
	libLapack.a(ieeeck.o)\
	libLapack.a(ilaenv.o)\
	libLapack.a(iparmq.o)\
	libLapack.a(zgebak.o)\
	libLapack.a(zgebal.o)\
	libLapack.a(zgeev.o)\
	libLapack.a(zgehd2.o)\
	libLapack.a(zgehrd.o)\
	libLapack.a(zgetf2.o)\
	libLapack.a(zgetrf.o)\
	libLapack.a(zhseqr.o)\
	libLapack.a(zlacgv.o)\
	libLapack.a(zlacpy.o)\
	libLapack.a(zladiv.o)\
	libLapack.a(zlahqr.o)\
	libLapack.a(zlahr2.o)\
	libLapack.a(zlange.o)\
	libLapack.a(zlaqr0.o)\
	libLapack.a(zlaqr1.o)\
	libLapack.a(zlaqr2.o)\
	libLapack.a(zlaqr3.o)\
	libLapack.a(zlaqr4.o)\
	libLapack.a(zlaqr5.o)\
	libLapack.a(zlarfb.o)\
	libLapack.a(zlarf.o)\
	libLapack.a(zlarfg.o)\
	libLapack.a(zlarft.o)\
	libLapack.a(zlartg.o)\
	libLapack.a(zlascl.o)\
	libLapack.a(zlaset.o)\
	libLapack.a(zlassq.o)\
	libLapack.a(zlaswp.o)\
	libLapack.a(zlatrs.o)\
	libLapack.a(zrot.o)\
	libLapack.a(ztrevc.o)\
	libLapack.a(ztrexc.o)\
	libLapack.a(zung2r.o)\
	libLapack.a(zunghr.o)\
	libLapack.a(zungqr.o) \
	libLapack.a(disnan.o )\
	libLapack.a(dlae2.o )\
	libLapack.a(dlaev2.o )\
	libLapack.a(dlaisnan.o )\
	libLapack.a(dlanst.o )\
	libLapack.a(dlansy.o )\
	libLapack.a(dlapy2.o )\
	libLapack.a(dlarf.o )\
	libLapack.a(dlarfb.o )\
	libLapack.a(dlarfg.o )\
	libLapack.a(dlarft.o )\
	libLapack.a(dlartg.o )\
	libLapack.a(dlascl.o )\
	libLapack.a(dlaset.o )\
	libLapack.a(dlasr.o )\
	libLapack.a(dlasrt.o )\
	libLapack.a(dlassq.o )\
	libLapack.a(dlatrd.o )\
	libLapack.a(dorg2l.o )\
	libLapack.a(dorg2r.o )\
	libLapack.a(dorgql.o )\
	libLapack.a(dorgqr.o )\
	libLapack.a(dorgtr.o )\
	libLapack.a(dsteqr.o )\
	libLapack.a(dsterf.o )\
	libLapack.a(dsyev.o )\
	libLapack.a(dsytd2.o )\
	libLapack.a(dsytrd.o )\
	libLapack.a(iladlr.o)\
	libLapack.a(iladlc.o)

