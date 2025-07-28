        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:22 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSBMV__genmod
          INTERFACE 
            SUBROUTINE DSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE DSBMV
          END INTERFACE 
        END MODULE DSBMV__genmod
