        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:12 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZHSEQR__genmod
          INTERFACE 
            SUBROUTINE ZHSEQR(JOB,COMPZ,N,ILO,IHI,H,LDH,W,Z,LDZ,WORK,   &
     &LWORK,INFO)
              INTEGER(KIND=4) :: LDZ
              INTEGER(KIND=4) :: LDH
              CHARACTER(LEN=1) :: JOB
              CHARACTER(LEN=1) :: COMPZ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              COMPLEX(KIND=8) :: H(LDH,*)
              COMPLEX(KIND=8) :: W(*)
              COMPLEX(KIND=8) :: Z(LDZ,*)
              COMPLEX(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZHSEQR
          END INTERFACE 
        END MODULE ZHSEQR__genmod
