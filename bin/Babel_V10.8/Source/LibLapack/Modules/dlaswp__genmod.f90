        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:11 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLASWP__genmod
          INTERFACE 
            SUBROUTINE DLASWP(N,A,LDA,K1,K2,IPIV,INCX)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: K1
              INTEGER(KIND=4) :: K2
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE DLASWP
          END INTERFACE 
        END MODULE DLASWP__genmod
