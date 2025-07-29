        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 20 15:48:17 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSYEV__genmod
          INTERFACE 
            SUBROUTINE DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: JOBZ
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: W(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSYEV
          END INTERFACE 
        END MODULE DSYEV__genmod
