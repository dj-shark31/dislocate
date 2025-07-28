        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:28 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZHBMV__genmod
          INTERFACE 
            SUBROUTINE ZHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE ZHBMV
          END INTERFACE 
        END MODULE ZHBMV__genmod
