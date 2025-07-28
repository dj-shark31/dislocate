        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:26 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STPSV__genmod
          INTERFACE 
            SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: AP(*)
              REAL(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE STPSV
          END INTERFACE 
        END MODULE STPSV__genmod
