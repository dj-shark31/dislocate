        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:24 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SDOT__genmod
          INTERFACE 
            FUNCTION SDOT(N,SX,INCX,SY,INCY)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: SX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: SY(*)
              INTEGER(KIND=4) :: INCY
              REAL(KIND=4) :: SDOT
            END FUNCTION SDOT
          END INTERFACE 
        END MODULE SDOT__genmod
