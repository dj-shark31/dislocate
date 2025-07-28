        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:21 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAXPY__genmod
          INTERFACE 
            SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: DA
              REAL(KIND=8) :: DX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: DY(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE DAXPY
          END INTERFACE 
        END MODULE DAXPY__genmod
