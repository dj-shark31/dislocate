        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:20 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CTPMV__genmod
          INTERFACE 
            SUBROUTINE CTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: AP(*)
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE CTPMV
          END INTERFACE 
        END MODULE CTPMV__genmod
