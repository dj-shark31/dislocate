        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:17 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ILADLC__genmod
          INTERFACE 
            FUNCTION ILADLC(M,N,A,LDA)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: ILADLC
            END FUNCTION ILADLC
          END INTERFACE 
        END MODULE ILADLC__genmod
