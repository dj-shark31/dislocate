        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:25 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SSPMV__genmod
          INTERFACE 
            SUBROUTINE SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: AP(*)
              REAL(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: BETA
              REAL(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE SSPMV
          END INTERFACE 
        END MODULE SSPMV__genmod
