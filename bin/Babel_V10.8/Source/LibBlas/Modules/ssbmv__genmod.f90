        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:25 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SSBMV__genmod
          INTERFACE 
            SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: A(LDA,*)
              REAL(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: BETA
              REAL(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE SSBMV
          END INTERFACE 
        END MODULE SSBMV__genmod
