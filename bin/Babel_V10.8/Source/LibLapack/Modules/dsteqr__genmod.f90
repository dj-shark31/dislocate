        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:17 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSTEQR__genmod
          INTERFACE 
            SUBROUTINE DSTEQR(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
              INTEGER(KIND=4) :: LDZ
              CHARACTER(LEN=1) :: COMPZ
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: Z(LDZ,*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSTEQR
          END INTERFACE 
        END MODULE DSTEQR__genmod
