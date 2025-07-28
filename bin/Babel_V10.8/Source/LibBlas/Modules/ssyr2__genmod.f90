        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:26 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SSYR2__genmod
          INTERFACE 
            SUBROUTINE SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
              REAL(KIND=4) :: A(LDA,*)
            END SUBROUTINE SSYR2
          END INTERFACE 
        END MODULE SSYR2__genmod
