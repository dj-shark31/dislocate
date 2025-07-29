        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:14 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZLASSQ__genmod
          INTERFACE 
            SUBROUTINE ZLASSQ(N,X,INCX,SCALE,SUMSQ)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: SCALE
              REAL(KIND=8) :: SUMSQ
            END SUBROUTINE ZLASSQ
          END INTERFACE 
        END MODULE ZLASSQ__genmod
