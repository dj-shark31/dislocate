        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:29 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZSYMM__genmod
          INTERFACE 
            SUBROUTINE ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: B(LDB,*)
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: C(LDC,*)
            END SUBROUTINE ZSYMM
          END INTERFACE 
        END MODULE ZSYMM__genmod
