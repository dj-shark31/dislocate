        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:18 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CGBMV__genmod
          INTERFACE 
            SUBROUTINE CGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y, &
     &INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: KL
              INTEGER(KIND=4) :: KU
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: A(LDA,*)
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=4) :: BETA
              COMPLEX(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE CGBMV
          END INTERFACE 
        END MODULE CGBMV__genmod
