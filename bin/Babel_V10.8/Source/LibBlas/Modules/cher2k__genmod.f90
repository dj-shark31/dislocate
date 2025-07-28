        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:19 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHER2K__genmod
          INTERFACE 
            SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,  &
     &LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: A(LDA,*)
              COMPLEX(KIND=4) :: B(LDB,*)
              REAL(KIND=4) :: BETA
              COMPLEX(KIND=4) :: C(LDC,*)
            END SUBROUTINE CHER2K
          END INTERFACE 
        END MODULE CHER2K__genmod
