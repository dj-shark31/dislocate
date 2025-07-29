        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:29 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZHPR__genmod
          INTERFACE 
            SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: AP(*)
            END SUBROUTINE ZHPR
          END INTERFACE 
        END MODULE ZHPR__genmod
