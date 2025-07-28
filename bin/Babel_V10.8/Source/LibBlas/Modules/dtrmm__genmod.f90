        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:23 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DTRMM__genmod
          INTERFACE 
            SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB&
     &)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANSA
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: B(LDB,*)
            END SUBROUTINE DTRMM
          END INTERFACE 
        END MODULE DTRMM__genmod
