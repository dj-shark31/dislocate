        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:25 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SSYMM__genmod
          INTERFACE 
            SUBROUTINE SSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: A(LDA,*)
              REAL(KIND=4) :: B(LDB,*)
              REAL(KIND=4) :: BETA
              REAL(KIND=4) :: C(LDC,*)
            END SUBROUTINE SSYMM
          END INTERFACE 
        END MODULE SSYMM__genmod
