        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:27 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STRSM__genmod
          INTERFACE 
            SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB&
     &)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANSA
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: A(LDA,*)
              REAL(KIND=4) :: B(LDB,*)
            END SUBROUTINE STRSM
          END INTERFACE 
        END MODULE STRSM__genmod
