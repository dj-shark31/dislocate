        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:22 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSDOT__genmod
          INTERFACE 
            FUNCTION DSDOT(N,SX,INCX,SY,INCY)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: SX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: SY(*)
              INTEGER(KIND=4) :: INCY
              REAL(KIND=8) :: DSDOT
            END FUNCTION DSDOT
          END INTERFACE 
        END MODULE DSDOT__genmod
