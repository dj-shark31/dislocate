        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:14 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZLASWP__genmod
          INTERFACE 
            SUBROUTINE ZLASWP(N,A,LDA,K1,K2,IPIV,INCX)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: K1
              INTEGER(KIND=4) :: K2
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE ZLASWP
          END INTERFACE 
        END MODULE ZLASWP__genmod
