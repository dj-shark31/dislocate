        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:28 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZHER2__genmod
          INTERFACE 
            SUBROUTINE ZHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
              COMPLEX(KIND=8) :: A(LDA,*)
            END SUBROUTINE ZHER2
          END INTERFACE 
        END MODULE ZHER2__genmod
