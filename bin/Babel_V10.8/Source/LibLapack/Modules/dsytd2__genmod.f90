        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:17 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSYTD2__genmod
          INTERFACE 
            SUBROUTINE DSYTD2(UPLO,N,A,LDA,D,E,TAU,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: D(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: TAU(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSYTD2
          END INTERFACE 
        END MODULE DSYTD2__genmod
