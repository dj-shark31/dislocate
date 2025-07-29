        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:23 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DTPMV__genmod
          INTERFACE 
            SUBROUTINE DTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: AP(*)
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE DTPMV
          END INTERFACE 
        END MODULE DTPMV__genmod
