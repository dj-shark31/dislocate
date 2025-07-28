        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:29 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZSYRK__genmod
          INTERFACE 
            SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: C(LDC,*)
            END SUBROUTINE ZSYRK
          END INTERFACE 
        END MODULE ZSYRK__genmod
