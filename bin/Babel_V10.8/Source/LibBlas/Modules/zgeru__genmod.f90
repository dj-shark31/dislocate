        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:28 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZGERU__genmod
          INTERFACE 
            SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
              COMPLEX(KIND=8) :: A(LDA,*)
            END SUBROUTINE ZGERU
          END INTERFACE 
        END MODULE ZGERU__genmod
