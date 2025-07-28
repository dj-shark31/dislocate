        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:11 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DGETF2__genmod
          INTERFACE 
            SUBROUTINE DGETF2(M,N,A,LDA,IPIV,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGETF2
          END INTERFACE 
        END MODULE DGETF2__genmod
