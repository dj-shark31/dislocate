        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:25 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SSPR__genmod
          INTERFACE 
            SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=4) :: AP(*)
            END SUBROUTINE SSPR
          END INTERFACE 
        END MODULE SSPR__genmod
