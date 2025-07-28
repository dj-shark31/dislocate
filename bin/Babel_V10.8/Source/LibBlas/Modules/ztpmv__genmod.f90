        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:30 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZTPMV__genmod
          INTERFACE 
            SUBROUTINE ZTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: AP(*)
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE ZTPMV
          END INTERFACE 
        END MODULE ZTPMV__genmod
