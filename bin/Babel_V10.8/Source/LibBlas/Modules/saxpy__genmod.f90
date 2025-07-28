        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:24 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SAXPY__genmod
          INTERFACE 
            SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: SA
              REAL(KIND=4) :: SX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: SY(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE SAXPY
          END INTERFACE 
        END MODULE SAXPY__genmod
