        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:13 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZLAQR5__genmod
          INTERFACE 
            SUBROUTINE ZLAQR5(WANTT,WANTZ,KACC22,N,KTOP,KBOT,NSHFTS,S,H,&
     &LDH,ILOZ,IHIZ,Z,LDZ,V,LDV,U,LDU,NV,WV,LDWV,NH,WH,LDWH)
              INTEGER(KIND=4) :: LDWH
              INTEGER(KIND=4) :: LDWV
              INTEGER(KIND=4) :: LDU
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: LDZ
              INTEGER(KIND=4) :: LDH
              LOGICAL(KIND=4) :: WANTT
              LOGICAL(KIND=4) :: WANTZ
              INTEGER(KIND=4) :: KACC22
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: KTOP
              INTEGER(KIND=4) :: KBOT
              INTEGER(KIND=4) :: NSHFTS
              COMPLEX(KIND=8) :: S(*)
              COMPLEX(KIND=8) :: H(LDH,*)
              INTEGER(KIND=4) :: ILOZ
              INTEGER(KIND=4) :: IHIZ
              COMPLEX(KIND=8) :: Z(LDZ,*)
              COMPLEX(KIND=8) :: V(LDV,*)
              COMPLEX(KIND=8) :: U(LDU,*)
              INTEGER(KIND=4) :: NV
              COMPLEX(KIND=8) :: WV(LDWV,*)
              INTEGER(KIND=4) :: NH
              COMPLEX(KIND=8) :: WH(LDWH,*)
            END SUBROUTINE ZLAQR5
          END INTERFACE 
        END MODULE ZLAQR5__genmod
