        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:29 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZSYR2K__genmod
          INTERFACE 
            SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,  &
     &LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: B(LDB,*)
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: C(LDC,*)
            END SUBROUTINE ZSYR2K
          END INTERFACE 
        END MODULE ZSYR2K__genmod
