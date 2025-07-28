        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:22 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSYMV__genmod
          INTERFACE 
            SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE DSYMV
          END INTERFACE 
        END MODULE DSYMV__genmod
