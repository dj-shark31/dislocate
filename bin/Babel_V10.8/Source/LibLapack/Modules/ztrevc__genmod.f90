        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:15 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZTREVC__genmod
          INTERFACE 
            SUBROUTINE ZTREVC(SIDE,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR&
     &,MM,M,WORK,RWORK,INFO)
              INTEGER(KIND=4) :: LDVR
              INTEGER(KIND=4) :: LDVL
              INTEGER(KIND=4) :: LDT
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: HOWMNY
              LOGICAL(KIND=4) :: SELECT(*)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: T(LDT,*)
              COMPLEX(KIND=8) :: VL(LDVL,*)
              COMPLEX(KIND=8) :: VR(LDVR,*)
              INTEGER(KIND=4) :: MM
              INTEGER(KIND=4) :: M
              COMPLEX(KIND=8) :: WORK(*)
              REAL(KIND=8) :: RWORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZTREVC
          END INTERFACE 
        END MODULE ZTREVC__genmod
