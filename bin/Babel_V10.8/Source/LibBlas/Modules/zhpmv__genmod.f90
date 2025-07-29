        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:29 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZHPMV__genmod
          INTERFACE 
            SUBROUTINE ZHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: AP(*)
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE ZHPMV
          END INTERFACE 
        END MODULE ZHPMV__genmod
