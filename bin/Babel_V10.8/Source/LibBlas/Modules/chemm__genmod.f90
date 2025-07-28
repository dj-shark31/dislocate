        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:18 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHEMM__genmod
          INTERFACE 
            SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: ALPHA
              COMPLEX(KIND=4) :: A(LDA,*)
              COMPLEX(KIND=4) :: B(LDB,*)
              COMPLEX(KIND=4) :: BETA
              COMPLEX(KIND=4) :: C(LDC,*)
            END SUBROUTINE CHEMM
          END INTERFACE 
        END MODULE CHEMM__genmod
