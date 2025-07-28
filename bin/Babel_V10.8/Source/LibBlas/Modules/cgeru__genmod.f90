        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:18 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CGERU__genmod
          INTERFACE 
            SUBROUTINE CGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=4) :: Y(*)
              INTEGER(KIND=4) :: INCY
              COMPLEX(KIND=4) :: A(LDA,*)
            END SUBROUTINE CGERU
          END INTERFACE 
        END MODULE CGERU__genmod
