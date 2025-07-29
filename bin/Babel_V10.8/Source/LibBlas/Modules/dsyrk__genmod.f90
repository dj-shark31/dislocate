        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:23 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSYRK__genmod
          INTERFACE 
            SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: C(LDC,*)
            END SUBROUTINE DSYRK
          END INTERFACE 
        END MODULE DSYRK__genmod
