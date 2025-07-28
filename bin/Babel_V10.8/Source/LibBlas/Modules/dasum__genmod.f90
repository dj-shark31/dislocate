        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:21 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DASUM__genmod
          INTERFACE 
            FUNCTION DASUM(N,DX,INCX)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: DX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: DASUM
            END FUNCTION DASUM
          END INTERFACE 
        END MODULE DASUM__genmod
