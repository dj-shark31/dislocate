        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:22 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSYR2K__genmod
          INTERFACE 
            SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,  &
     &LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: C(LDC,*)
            END SUBROUTINE DSYR2K
          END INTERFACE 
        END MODULE DSYR2K__genmod
