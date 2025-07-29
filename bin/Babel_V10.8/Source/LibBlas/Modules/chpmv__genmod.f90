        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:19 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHPMV__genmod
          INTERFACE 
            SUBROUTINE CHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: AP(*)
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=4) :: BETA
              COMPLEX(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE CHPMV
          END INTERFACE 
        END MODULE CHPMV__genmod
