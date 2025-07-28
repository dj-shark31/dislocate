.PRECIOUS: libBlas.a
.SUFFIXES: .f .o .a
MODULES = ./Modules


.f.o:
	$(F77) $(F77FLAGS) -c $<

.f.a:
	$(F77) $(F77FLAGS) -c $<
	ar Urv $@ $%
	@rm -f $%

all:	libBlas.a

clean:
	@rm -f *.o *.mod $(MODULES)/*.mod libBlas.a

libBlas.a: libBlas.a(caxpy.o) \
	libBlas.a(ccopy.o) \
	libBlas.a(cdotc.o) \
	libBlas.a(cdotu.o) \
	libBlas.a(cgbmv.o) \
	libBlas.a(cgemm.o) \
	libBlas.a(cgemv.o) \
	libBlas.a(cgerc.o) \
	libBlas.a(cgeru.o) \
	libBlas.a(chbmv.o) \
	libBlas.a(chemm.o) \
	libBlas.a(chemv.o) \
	libBlas.a(cher2.o) \
	libBlas.a(cher2k.o) \
	libBlas.a(cher.o) \
	libBlas.a(cherk.o) \
	libBlas.a(chpmv.o) \
	libBlas.a(chpr2.o) \
	libBlas.a(chpr.o) \
	libBlas.a(crotg.o) \
	libBlas.a(cscal.o) \
	libBlas.a(csrot.o) \
	libBlas.a(csscal.o) \
	libBlas.a(cswap.o) \
	libBlas.a(csymm.o) \
	libBlas.a(csyr2k.o) \
	libBlas.a(csyrk.o) \
	libBlas.a(ctbmv.o) \
	libBlas.a(ctbsv.o) \
	libBlas.a(ctpmv.o) \
	libBlas.a(ctpsv.o) \
	libBlas.a(ctrmm.o) \
	libBlas.a(ctrmv.o) \
	libBlas.a(ctrsm.o) \
	libBlas.a(ctrsv.o) \
	libBlas.a(dasum.o) \
	libBlas.a(daxpy.o) \
	libBlas.a(dcabs1.o) \
	libBlas.a(dcopy.o) \
	libBlas.a(ddot.o) \
	libBlas.a(dgbmv.o) \
	libBlas.a(dgemm.o) \
	libBlas.a(dgemv.o) \
	libBlas.a(dger.o) \
	libBlas.a(dnrm2.o) \
	libBlas.a(drot.o) \
	libBlas.a(drotg.o) \
	libBlas.a(drotm.o) \
	libBlas.a(drotmg.o) \
	libBlas.a(dsbmv.o) \
	libBlas.a(dscal.o) \
	libBlas.a(dsdot.o) \
	libBlas.a(dspmv.o) \
	libBlas.a(dspr2.o) \
	libBlas.a(dspr.o) \
	libBlas.a(dswap.o) \
	libBlas.a(dsymm.o) \
	libBlas.a(dsymv.o) \
	libBlas.a(dsyr2.o) \
	libBlas.a(dsyr2k.o) \
	libBlas.a(dsyr.o) \
	libBlas.a(dsyrk.o) \
	libBlas.a(dtbmv.o) \
	libBlas.a(dtbsv.o) \
	libBlas.a(dtpmv.o) \
	libBlas.a(dtpsv.o) \
	libBlas.a(dtrmm.o) \
	libBlas.a(dtrmv.o) \
	libBlas.a(dtrsm.o) \
	libBlas.a(dtrsv.o) \
	libBlas.a(dzasum.o) \
	libBlas.a(dznrm2.o) \
	libBlas.a(icamax.o) \
	libBlas.a(idamax.o) \
	libBlas.a(isamax.o) \
	libBlas.a(izamax.o) \
	libBlas.a(lsame.o) \
	libBlas.a(sasum.o) \
	libBlas.a(saxpy.o) \
	libBlas.a(scasum.o) \
	libBlas.a(scnrm2.o) \
	libBlas.a(scopy.o) \
	libBlas.a(sdot.o) \
	libBlas.a(sdsdot.o) \
	libBlas.a(sgbmv.o) \
	libBlas.a(sgemm.o) \
	libBlas.a(sgemv.o) \
	libBlas.a(sger.o) \
	libBlas.a(snrm2.o) \
	libBlas.a(srot.o) \
	libBlas.a(srotg.o) \
	libBlas.a(srotm.o) \
	libBlas.a(srotmg.o) \
	libBlas.a(ssbmv.o) \
	libBlas.a(sscal.o) \
	libBlas.a(sspmv.o) \
	libBlas.a(sspr2.o) \
	libBlas.a(sspr.o) \
	libBlas.a(sswap.o) \
	libBlas.a(ssymm.o) \
	libBlas.a(ssymv.o) \
	libBlas.a(ssyr2.o) \
	libBlas.a(ssyr2k.o) \
	libBlas.a(ssyr.o) \
	libBlas.a(ssyrk.o) \
	libBlas.a(stbmv.o) \
	libBlas.a(stbsv.o) \
	libBlas.a(stpmv.o) \
	libBlas.a(stpsv.o) \
	libBlas.a(strmm.o) \
	libBlas.a(strmv.o) \
	libBlas.a(strsm.o) \
	libBlas.a(strsv.o) \
	libBlas.a(xerbla.o) \
	libBlas.a(zaxpy.o) \
	libBlas.a(zcopy.o) \
	libBlas.a(zdotc.o) \
	libBlas.a(zdotu.o) \
	libBlas.a(zdrot.o) \
	libBlas.a(zdscal.o) \
	libBlas.a(zgbmv.o) \
	libBlas.a(zgemm.o) \
	libBlas.a(zgemv.o) \
	libBlas.a(zgerc.o) \
	libBlas.a(zgeru.o) \
	libBlas.a(zhbmv.o) \
	libBlas.a(zhemm.o) \
	libBlas.a(zhemv.o) \
	libBlas.a(zher2.o) \
	libBlas.a(zher2k.o) \
	libBlas.a(zher.o) \
	libBlas.a(zherk.o) \
	libBlas.a(zhpmv.o) \
	libBlas.a(zhpr2.o) \
	libBlas.a(zhpr.o) \
	libBlas.a(zrotg.o) \
	libBlas.a(zscal.o) \
	libBlas.a(zswap.o) \
	libBlas.a(zsymm.o) \
	libBlas.a(zsyr2k.o) \
	libBlas.a(zsyrk.o) \
	libBlas.a(ztbmv.o) \
	libBlas.a(ztbsv.o) \
	libBlas.a(ztpmv.o) \
	libBlas.a(ztpsv.o) \
	libBlas.a(ztrmm.o) \
	libBlas.a(ztrmv.o) \
	libBlas.a(ztrsm.o) \
	libBlas.a(ztrsv.o)
