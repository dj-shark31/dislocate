        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:19 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHPR__genmod
          INTERFACE 
            SUBROUTINE CHPR(UPLO,N,ALPHA,X,INCX,AP)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=4) :: AP(*)
            END SUBROUTINE CHPR
          END INTERFACE 
        END MODULE CHPR__genmod
