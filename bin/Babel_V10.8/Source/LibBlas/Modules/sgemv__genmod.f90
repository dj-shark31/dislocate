        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:24 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SGEMV__genmod
          INTERFACE 
            SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: A(LDA,*)
              REAL(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: BETA
              REAL(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE SGEMV
          END INTERFACE 
        END MODULE SGEMV__genmod
