        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:13 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZLAQR2__genmod
          INTERFACE 
            SUBROUTINE ZLAQR2(WANTT,WANTZ,N,KTOP,KBOT,NW,H,LDH,ILOZ,IHIZ&
     &,Z,LDZ,NS,ND,SH,V,LDV,NH,T,LDT,NV,WV,LDWV,WORK,LWORK)
              INTEGER(KIND=4) :: LDWV
              INTEGER(KIND=4) :: LDT
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: LDZ
              INTEGER(KIND=4) :: LDH
              LOGICAL(KIND=4) :: WANTT
              LOGICAL(KIND=4) :: WANTZ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: KTOP
              INTEGER(KIND=4) :: KBOT
              INTEGER(KIND=4) :: NW
              COMPLEX(KIND=8) :: H(LDH,*)
              INTEGER(KIND=4) :: ILOZ
              INTEGER(KIND=4) :: IHIZ
              COMPLEX(KIND=8) :: Z(LDZ,*)
              INTEGER(KIND=4) :: NS
              INTEGER(KIND=4) :: ND
              COMPLEX(KIND=8) :: SH(*)
              COMPLEX(KIND=8) :: V(LDV,*)
              INTEGER(KIND=4) :: NH
              COMPLEX(KIND=8) :: T(LDT,*)
              INTEGER(KIND=4) :: NV
              COMPLEX(KIND=8) :: WV(LDWV,*)
              COMPLEX(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
            END SUBROUTINE ZLAQR2
          END INTERFACE 
        END MODULE ZLAQR2__genmod
