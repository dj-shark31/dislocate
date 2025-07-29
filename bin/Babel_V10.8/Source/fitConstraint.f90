MODULE FitConstraintModule

  USE disloc_elasticity_ani, ONLY : max_nd
  USE dDipoleModule, ONLY : max_nDDipole
  IMPLICIT NONE
  SAVE

  ! Variable set to true if constraint are applied
  LOGICAL :: l_constraint

  ! Number of constraints and maximum value
  INTEGER :: nCons
  INTEGER, parameter :: max_nCons=20

  ! Fixed value for each constraint
  REAL(kind(0.d0)), dimension(1:max_nCons), private :: consValue

  ! cons_bDislo(nc,i,n)=.true. if bDislo(i,n) takes part in constraint nc
  LOGICAL, dimension(1:max_nCons,1:3,1:max_nd), private :: cons_bDislo
  LOGICAL, dimension(1:max_nCons,1:3,1:max_nDDipole), private :: cons_bDDipole

  REAL(kind(0.d0)), dimension(1:max_nCons,1:3,1:max_nd), private :: weight_bDislo
  REAL(kind(0.d0)), dimension(1:max_nCons,1:3,1:max_nDDipole), private :: weight_bDDipole

  PRIVATE :: calculateConstraint

CONTAINS

  SUBROUTINE InitConstraint()

    IMPLICIT NONE

    l_constraint = .FALSE.
    nCons = 0
    consValue(:) = 0.d0
    cons_bDislo(:,:,:) = .FALSE.
    cons_bDDipole(:,:,:) = .FALSE.
    weight_bDislo(:,:,:) = 1.d0
    weight_bDDipole(:,:,:) = 1.d0

  END SUBROUTINE InitConstraint

  SUBROUTINE ReadConstraint(inp)
  
    IMPLICIT NONE
    INTEGER, intent(in) :: inp  ! Input unit

    INTEGER, dimension(1:3,1:max_nd) :: bDislo
    INTEGER, dimension(1:3,1:max_nDDipole) :: bDDipole
    INTEGER :: n, i
    REAL(kind(0.d0)) :: weight_sum

    NAMELIST /constraint/ bDislo, bDDipole

    ! Initialization
    bDislo(:,:) = 0
    bDDipole(:,:) = 0

    ! Read input data
    READ(inp,nml=constraint)

    ! Increment number of constraints
    nCons = nCons + 1
    l_constraint = .true.

    IF (nCons.GT.max_nCons) THEN
            WRITE(0,'(a)') "Value of parameter max_nCons too small"
            STOP "< ReadConstraint >"
    END IF

    ! bDislo(:,:)=0 => variable not constrained
    cons_bDislo(nCons,:,:) = bDislo(:,:).NE.0
    cons_bDDipole(nCons,:,:) = bDDipole(:,:).NE.0

    ! Calculate fixed value for the constraint
    consValue(nCons) = calculateConstraint(nCons)

    ! Check that constraint can be applied
    weight_sum = Sum( weight_bDislo(nCons,:,:), cons_bDislo(nCons,:,:) ) &
                + Sum( weight_bDDipole(nCons,:,:), cons_bDDipole(nCons,:,:) )
    IF (Abs(weight_sum).LE.1.d-12) THEN
            WRITE(0,'(a,i0,a)') 'Constraint ', nCons, ' cannot be applied as weight sum equals 0'
            CALL PrintConstraint(0, nCons) 
            STOP '< ReadConstraint >'
    END IF

    ! Modify weights associated with constrained Burgers vectors
    ! in order for the quantity not to be modified by other constraints
    DO n=1, max_nd
       DO i=1, 3
          IF (cons_bDislo(nCons,i,n)) weight_bDislo(nCons+1:max_nCons,i,n)=0.d0
       END DO
    END DO
    DO n=1, max_nDDipole
       DO i=1, 3
          IF (cons_bDDipole(nCons,i,n)) weight_bDDipole(nCons+1:max_nCons,i,n)=0.d0
       END DO
    END DO


  END SUBROUTINE ReadConstraint

  FUNCTION calculateConstraint(nc) RESULT(bSum)

    USE disloc_elasticity_ani, ONLY : nd, bDislo
    USE dDipoleModule, ONLY : nDDipole, bDDipole
    IMPLICIT NONE
    INTEGER, intent(in) :: nc
    REAL(kind(0.d0)) :: bSum

    bSum = Sum( bDislo(:,1:nd), cons_bDislo(nc,:,1:nd) ) &
        + Sum( bDDipole(:,1:nDDipole), cons_bDDipole(nc,:,1:nDDipole) )

  END FUNCTION calculateConstraint

  SUBROUTINE ScaleConstraint(a)
    ! Multiply distances by a for constraint definition

    IMPLICIT NONE
    REAL(kind(0.d0)), intent(in) :: a
    INTEGER :: nc

    DO nc=1, nCons
      consValue(nc) = a*consValue(nc)
    END DO

  END SUBROUTINE ScaleConstraint
  
  SUBROUTINE PrintConstraint(out, ncOpt)
  
    USE disloc_elasticity_ani, ONLY : nd
    USE dDipoleModule, ONLY : nDDipole
    IMPLICIT NONE

    INTEGER, intent(in) :: out
    INTEGER, intent(in), optional :: ncOpt      ! Only this constraint is printed

    CHARACTER(len=15) :: phrase
    INTEGER :: nc, n

    WRITE(out,*)
    WRITE(out,'(a,i0)') 'Fit: number of constraints: ', nCons

    DO nc=1, nCons
       IF (present(ncOpt)) THEN
               IF (nc.NE.ncOpt) Cycle
       END IF
       WRITE(out,'(a,i0)') '  Constraint ', nc
       DO n=1, nd
          IF (Count(cons_bDislo(nc,:,n)).GT.0) THEN
                  phrase=''
                  IF (cons_bDislo(nc,1,n)) WRITE(phrase,'(2a)') Trim(phrase), ' bx'
                  IF (cons_bDislo(nc,2,n)) WRITE(phrase,'(2a)') Trim(phrase), ' by'
                  IF (cons_bDislo(nc,3,n)) WRITE(phrase,'(2a)') Trim(phrase), ' bz'
                  WRITE(out,'(a,i0,3a)') '    dislo ', n, ':', phrase, ' constrained'
          END IF
       END DO
       DO n=1, nDDipole
          IF (Count(cons_bDDipole(nc,:,n)).GT.0) THEN
                  phrase=''
                  IF (cons_bDDipole(nc,1,n)) WRITE(phrase,'(2a)') Trim(phrase), ' bx'
                  IF (cons_bDDipole(nc,2,n)) WRITE(phrase,'(2a)') Trim(phrase), ' by'
                  IF (cons_bDDipole(nc,3,n)) WRITE(phrase,'(2a)') Trim(phrase), ' bz'
                  WRITE(out,'(a,i0,3a)') '    dipole ', n, ':', phrase, ' constrained'
          END IF
       END DO
       WRITE(out,'(a,g20.12)') '    Sum should equal ', consValue(nc)
    END DO
    WRITE(out,*)

  END SUBROUTINE PrintConstraint

  SUBROUTINE ApplyConstraint()

    USE disloc_elasticity_ani, ONLY : nd, bDislo
    USE dDipoleModule, ONLY : nDDipole, bDDipole
    USE rearrange
    IMPLICIT NONE

    INTEGER :: nc, n, n2, i
    REAL(kind(0.d0)) :: bSum, bCorr

    ! Loop on constraints
    DO nc=1, nCons

      ! Calculate actual sum
      bSum = Sum( bDislo(:,1:nd), cons_bDislo(nc,:,1:nd) ) &
        + Sum( bDDipole(:,1:nDDipole), cons_bDDipole(nc,:,1:nDDipole) )

      ! Correction to apply to each Burgers vector
      !!bCorr = ( consValue(nc) - bSum )/Dble( Count( cons_bDislo(nc,:,1:nd) ) &
                !!+ Count( cons_bDDipole(nc,:,1:nDDipole) ) )
      bCorr = ( consValue(nc) - bSum )/( Sum( weight_bDislo(nc,:,1:nd), cons_bDislo(nc,:,1:nd) ) &
                + Sum( weight_bDDipole(nc,:,1:nDDipole), cons_bDDipole(nc,:,1:nDDipole) ) )


      ! Apply correction to Burgers vector to obtain the correct sum
      DO n=1, nd
         DO i=1, 3
            IF (cons_bDislo(nc,i,n)) THEN
                    !!bDislo(i,n) = bDislo(i,n) + bCorr
                    bDislo(i,n) = bDislo(i,n) + weight_bDislo(nc,i,n)*bCorr
                    ! Initialize Burgers vector of second dislocation composing the dipole
                    n2 = iDisloDipole(n)
                    IF (n2.GT.0) bDislo(i,n2) = -bDislo(i,n)
            END IF
         END DO
      END DO
      DO n=1, nDDipole
         DO i=1, 3
            IF (cons_bDDipole(nc,i,n)) THEN
                    !!bDDipole(i,n) = bDDipole(i,n) + bCorr
                    bDDipole(i,n) = bDDipole(i,n) + weight_bDDipole(nc,i,n)*bCorr
            END IF
         END DO
      END DO

    END DO

  END SUBROUTINE ApplyConstraint

END MODULE FitConstraintModule
