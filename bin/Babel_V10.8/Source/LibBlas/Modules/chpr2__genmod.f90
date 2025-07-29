        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:19 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHPR2__genmod
          INTERFACE 
            SUBROUTINE CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
              COMPLEX(KIND=4) :: AP(*)
            END SUBROUTINE CHPR2
          END INTERFACE 
        END MODULE CHPR2__genmod
