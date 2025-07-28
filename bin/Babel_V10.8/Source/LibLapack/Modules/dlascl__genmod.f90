        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:16 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLASCL__genmod
          INTERFACE 
            SUBROUTINE DLASCL(TYPE,KL,KU,CFROM,CTO,M,N,A,LDA,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TYPE
              INTEGER(KIND=4) :: KL
              INTEGER(KIND=4) :: KU
              REAL(KIND=8) :: CFROM
              REAL(KIND=8) :: CTO
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DLASCL
          END INTERFACE 
        END MODULE DLASCL__genmod
