        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:20 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CTRSV__genmod
          INTERFACE 
            SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: A(LDA,*)
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE CTRSV
          END INTERFACE 
        END MODULE CTRSV__genmod
