        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:20 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CSYRK__genmod
          INTERFACE 
            SUBROUTINE CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: A(LDA,*)
              COMPLEX(KIND=4) :: BETA
              COMPLEX(KIND=4) :: C(LDC,*)
            END SUBROUTINE CSYRK
          END INTERFACE 
        END MODULE CSYRK__genmod
