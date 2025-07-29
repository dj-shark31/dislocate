        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:12 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZLAHQR__genmod
          INTERFACE 
            SUBROUTINE ZLAHQR(WANTT,WANTZ,N,ILO,IHI,H,LDH,W,ILOZ,IHIZ,Z,&
     &LDZ,INFO)
              INTEGER(KIND=4) :: LDZ
              INTEGER(KIND=4) :: LDH
              LOGICAL(KIND=4) :: WANTT
              LOGICAL(KIND=4) :: WANTZ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              COMPLEX(KIND=8) :: H(LDH,*)
              COMPLEX(KIND=8) :: W(*)
              INTEGER(KIND=4) :: ILOZ
              INTEGER(KIND=4) :: IHIZ
              COMPLEX(KIND=8) :: Z(LDZ,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZLAHQR
          END INTERFACE 
        END MODULE ZLAHQR__genmod
