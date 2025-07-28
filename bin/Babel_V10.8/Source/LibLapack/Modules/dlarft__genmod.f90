        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:16 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLARFT__genmod
          INTERFACE 
            SUBROUTINE DLARFT(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
              INTEGER(KIND=4) :: LDT
              INTEGER(KIND=4) :: LDV
              CHARACTER(LEN=1) :: DIRECT
              CHARACTER(LEN=1) :: STOREV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: T(LDT,*)
            END SUBROUTINE DLARFT
          END INTERFACE 
        END MODULE DLARFT__genmod
