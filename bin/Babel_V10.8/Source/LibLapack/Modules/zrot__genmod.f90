        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:15 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZROT__genmod
          INTERFACE 
            SUBROUTINE ZROT(N,CX,INCX,CY,INCY,C,S)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: CX(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: CY(*)
              INTEGER(KIND=4) :: INCY
              REAL(KIND=8) :: C
              COMPLEX(KIND=8) :: S
            END SUBROUTINE ZROT
          END INTERFACE 
        END MODULE ZROT__genmod
