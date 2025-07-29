        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:12 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZGEEV__genmod
          INTERFACE 
            SUBROUTINE ZGEEV(JOBVL,JOBVR,N,A,LDA,W,VL,LDVL,VR,LDVR,WORK,&
     &LWORK,RWORK,INFO)
              INTEGER(KIND=4) :: LDVR
              INTEGER(KIND=4) :: LDVL
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: JOBVL
              CHARACTER(LEN=1) :: JOBVR
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: W(*)
              COMPLEX(KIND=8) :: VL(LDVL,*)
              COMPLEX(KIND=8) :: VR(LDVR,*)
              COMPLEX(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              REAL(KIND=8) :: RWORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZGEEV
          END INTERFACE 
        END MODULE ZGEEV__genmod
