        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:16 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLASR__genmod
          INTERFACE 
            SUBROUTINE DLASR(SIDE,PIVOT,DIRECT,M,N,C,S,A,LDA)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: PIVOT
              CHARACTER(LEN=1) :: DIRECT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: A(LDA,*)
            END SUBROUTINE DLASR
          END INTERFACE 
        END MODULE DLASR__genmod
