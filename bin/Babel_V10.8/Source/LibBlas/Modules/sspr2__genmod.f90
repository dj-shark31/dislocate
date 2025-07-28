        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:25 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SSPR2__genmod
          INTERFACE 
            SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
              REAL(KIND=4) :: AP(*)
            END SUBROUTINE SSPR2
          END INTERFACE 
        END MODULE SSPR2__genmod
