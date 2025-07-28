        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:22 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSPR2__genmod
          INTERFACE 
            SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
              REAL(KIND=8) :: AP(*)
            END SUBROUTINE DSPR2
          END INTERFACE 
        END MODULE DSPR2__genmod
