        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:11 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZGEBAK__genmod
          INTERFACE 
            SUBROUTINE ZGEBAK(JOB,SIDE,N,ILO,IHI,SCALE,M,V,LDV,INFO)
              INTEGER(KIND=4) :: LDV
              CHARACTER(LEN=1) :: JOB
              CHARACTER(LEN=1) :: SIDE
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: SCALE(*)
              INTEGER(KIND=4) :: M
              COMPLEX(KIND=8) :: V(LDV,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZGEBAK
          END INTERFACE 
        END MODULE ZGEBAK__genmod
