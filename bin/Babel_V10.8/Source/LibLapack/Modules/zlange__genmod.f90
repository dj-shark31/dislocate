        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:13 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZLANGE__genmod
          INTERFACE 
            FUNCTION ZLANGE(NORM,M,N,A,LDA,WORK)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: NORM
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: WORK(*)
              REAL(KIND=8) :: ZLANGE
            END FUNCTION ZLANGE
          END INTERFACE 
        END MODULE ZLANGE__genmod
