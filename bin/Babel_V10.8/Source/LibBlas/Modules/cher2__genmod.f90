        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:19 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHER2__genmod
          INTERFACE 
            SUBROUTINE CHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
              COMPLEX(KIND=4) :: A(LDA,*)
            END SUBROUTINE CHER2
          END INTERFACE 
        END MODULE CHER2__genmod
