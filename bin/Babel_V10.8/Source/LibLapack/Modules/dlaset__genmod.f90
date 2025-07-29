        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:16 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLASET__genmod
          INTERFACE 
            SUBROUTINE DLASET(UPLO,M,N,ALPHA,BETA,A,LDA)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: A(LDA,*)
            END SUBROUTINE DLASET
          END INTERFACE 
        END MODULE DLASET__genmod
