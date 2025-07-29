        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:24 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SCASUM__genmod
          INTERFACE 
            FUNCTION SCASUM(N,CX,INCX)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: CX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: SCASUM
            END FUNCTION SCASUM
          END INTERFACE 
        END MODULE SCASUM__genmod
