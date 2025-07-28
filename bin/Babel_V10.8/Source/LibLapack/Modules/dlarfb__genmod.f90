        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:16 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLARFB__genmod
          INTERFACE 
            SUBROUTINE DLARFB(SIDE,TRANS,DIRECT,STOREV,M,N,K,V,LDV,T,LDT&
     &,C,LDC,WORK,LDWORK)
              INTEGER(KIND=4) :: LDWORK
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDT
              INTEGER(KIND=4) :: LDV
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIRECT
              CHARACTER(LEN=1) :: STOREV
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: T(LDT,*)
              REAL(KIND=8) :: C(LDC,*)
              REAL(KIND=8) :: WORK(LDWORK,*)
            END SUBROUTINE DLARFB
          END INTERFACE 
        END MODULE DLARFB__genmod
