        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:24 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SGEMM__genmod
          INTERFACE 
            SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,&
     &C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANSA
              CHARACTER(LEN=1) :: TRANSB
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: A(LDA,*)
              REAL(KIND=4) :: B(LDB,*)
              REAL(KIND=4) :: BETA
              REAL(KIND=4) :: C(LDC,*)
            END SUBROUTINE SGEMM
          END INTERFACE 
        END MODULE SGEMM__genmod
