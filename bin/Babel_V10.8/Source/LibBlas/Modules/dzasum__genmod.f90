        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:23 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DZASUM__genmod
          INTERFACE 
            FUNCTION DZASUM(N,ZX,INCX)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ZX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: DZASUM
            END FUNCTION DZASUM
          END INTERFACE 
        END MODULE DZASUM__genmod
