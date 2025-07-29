        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:27 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZCOPY__genmod
          INTERFACE 
            SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ZX(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: ZY(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE ZCOPY
          END INTERFACE 
        END MODULE ZCOPY__genmod
