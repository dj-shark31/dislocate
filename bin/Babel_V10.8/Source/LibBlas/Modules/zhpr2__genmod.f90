        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:29 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZHPR2__genmod
          INTERFACE 
            SUBROUTINE ZHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
              COMPLEX(KIND=8) :: AP(*)
            END SUBROUTINE ZHPR2
          END INTERFACE 
        END MODULE ZHPR2__genmod
