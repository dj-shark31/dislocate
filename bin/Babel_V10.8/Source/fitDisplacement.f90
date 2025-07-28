MODULE fitDisplacementModule

  IMPLICIT NONE
  SAVE

  ! Number of atoms
  INTEGER, private :: imDis
  ! Coordinates of each atoms
  REAL(kind(0.d0)), dimension(:,:), allocatable, private :: xpDis
  ! Displacement of each atom
  REAL(kind(0.d0)), dimension(:,:), allocatable, private :: u0Dis, uDis
  ! Name of each atom
  CHARACTER(len=5), dimension(:), allocatable, private :: labelDis

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE ReadDisplacement(inp)
    ! Read displacement in file connected to unit inp
    ! Only points for which the differential displacement is greater than
    ! threshold are kept

    IMPLICIT NONE
    INTEGER, intent(in) :: inp  ! Unit connected to displacement file

    INTEGER :: i

    imDis=0

    
    READ(inp,*) imDis   ! Number of atoms
    READ(inp,*)         ! Line for title
    ALLOCATE(xpDis(1:3,1:imDis))
    ALLOCATE(u0Dis(1:3,1:imDis))
    ALLOCATE(labelDis(1:imDis))
    DO i=1, imDis
       READ(inp,*) labelDis(i), xpDis(1:3,i), u0Dis(1:3,i)
    END DO

  END SUBROUTINE ReadDisplacement

  SUBROUTINE WriteDisplacement(uLocal,out)

    IMPLICIT NONE
    INTEGER, intent(in) :: out  ! Unit connected to displacement file
    REAL(kind(0.d0)), dimension(:,:), intent(in) :: uLocal

    INTEGER :: i

    WRITE(out,'(i0)') imDis
    WRITE(out,*)
    DO i=1, imDis
       WRITE(out,'(a,6g20.12)') labelDis(i), xpDis(1:3,i), uLocal(1:3,i)
    END DO

  END SUBROUTINE WriteDisplacement

  SUBROUTINE FitDisplacement(out)

    USE Babel_data
    USE Math
    USE PeriodRelaxModule
    USE Powell_min
    USE FitModule
    USE FitConstraintModule
    IMPLICIT NONE
    INTEGER, intent(in) :: out

    REAL(kind(0.d0)), dimension(1:3,1:3) :: epsi
    INTEGER :: i, n, nDim, iter
    REAL(kind(0.d0)), dimension(:), allocatable :: x
    REAL(kind(0.d0)), dimension(:,:), allocatable :: xi
    REAL(kind(0.d0)) :: fCost

    ! ==== Fix homogeneous strain ========
    ! Homogeneous strain
    IF (.NOT.strain) THEN
            eStrain = 0.d0
    END IF

    ! Add homogeneous strain induced by line defects
    IF (induced_homogeneous_strain) THEN
            CALL RelaxPeriod(epsi,out)
            eStrain = eStrain + epsi
    END IF

    IF (matNorm2(eStrain).GT.1.d-6) strain=.true.
    ! Corresponding stress
    sStrainVoigt(1:6) = MatMul( CVoigt, (/ eStrain(1,1), eStrain(2,2), eStrain(3,3), &
        eStrain(2,3)+eStrain(3,2), eStrain(1,3)+eStrain(3,1), eStrain(1,2)+eStrain(2,1) /) )
    IF (verbosity.GE.verbosity_max) THEN
            IF (xImages.OR.yImages.OR.zImages) THEN
                    WRITE(out,*)
                    WRITE(out,'(a)') '--------------'
                    WRITE(out,*) 
                    WRITE(out,'(a)') 'Homogeneous strain, applied and resulting&
                        & from all line-defects'
                    DO i=1,3
                       WRITE(out,'(a,i1,a,3(g11.4,2x),a)') '  e(',i,',1:3) = | ', &
                            eStrain(i,1:3), ' |'
                    END DO
            END IF
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
    END IF

    ! ==== Initialize variable to dislocation positions
    nDim = nDimFit()
    Allocate(x(1:nDim))
    CALL Fit_WriteX(x, nDim)

    ! ==== Allocate array used by the cost function
    ALLOCATE(uDis(1:3,1:imDis))

    ! ==== First call =====
    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Minimization of Displacement'
            WRITE(out,'(a,i0)') '  nDim = ', nDim
            WRITE(out,*)
    END IF

    !==== Minimization =================
    Allocate(xi(1:nDim,1:nDim))
    xi(:,:)=0.d0
    DO n=1, nDim
       xi(n,n)=1.d0
    END DO
    iter = 0
    CALL powell(DisplacementCost,x,xi,nDim,nDim,ftol,iter,fCost, print_DisplacementCost)

    ! ==== Get result ==========
    CALL Fit_ReadX(x, nDim)
    ! ... and apply constraints
    CALL ApplyConstraint()

    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'End of minimization'
            CALL Fit_PrintX(x, nDim, out)
            fCost = DisplacementCost(x, nDim)
            WRITE(out,'(a,g20.12)') '  fcost = ', fcost
            WRITE(out,*)
    END IF

    DeAllocate(xi)
    DeAllocate(x)

    ! ==== Write obtained Displacement map in a file =======
    OPEN(file='displacement.fit',unit=60,action='write')
    CALL WriteDisplacement(uDis,60)
    CLOSE(60)
    OPEN(file='displacement.residual',unit=61,action='write')
    CALL WriteDisplacement(u0Dis-uDis,61)
    CLOSE(61)
    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Obtained displacement map written in file displacement.fit'
            WRITE(out,'(a)') 'Difference between input and fitted displacement map written in file displacement.residual'
            WRITE(out,*)
    END IF

    ! ==== The End =====
    DEALLOCATE(uDis)

  END SUBROUTINE FitDisplacement

  FUNCTION DisplacementCost(x, nDim) RESULT(fCost)

    USE Babel_data
    USE Math
    USE Disloc_elasticity_ani
    USE DDipoleModule
    USE PeriodModule
    USE rearrange
    USE FitModule
    USE FitConstraintModule
    IMPLICIT NONE
    ! Number of variables (3*number of dislocations)
    INTEGER, intent(in) :: nDim
    ! Variables (dislocation positions)
    REAL(kind(0.d0)), dimension(1:nDim), intent(in) :: x
    ! Cost function
    REAL(kind(0.d0)) :: fCost

    INTEGER :: i, n
    INTEGER :: out, verbosity_backup

    ! ==== Output unit and verbosity =======
    out=6
    IF (verbosity.LT.verbosity_max+1) THEN
            verbosity_backup=verbosity
            verbosity=0
    ELSE
            verbosity_backup=-1
    END IF

    IF (verbosity.GE.verbosity_max+1) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'DisplacementCost CALL'
            WRITE(out,'(a,i0)') '  nDim = ', nDim
            DO n=1, nDim
               WRITE(out,'(a,i0,a,g20.12)') '  x(',n,') = ', x(n)
            END DO
            WRITE(out,*)
    END IF


    ! ==== Initialize dislocation positions 
    CALL Fit_ReadX(x, nDim)
    ! ... and apply constraints
    CALL ApplyConstraint()

    ! ==== Initialize dislocation cuts and elastic calculation ========
    IF (l_dislo) THEN             ! Initialization for dislocation
            IF (verbosity.GE.verbosity_max+1) THEN
                    WRITE(out,'(a)') 'Initialization for dislocation'
                    WRITE(out,*)
            END IF
            CALL RearrangeDislo(out)
            CALL InitDisloc_ani(out)
    END IF
    ! ==== Initialize dislocation cuts and elastic calculation ========
    IF (l_DDipole) THEN             ! Initialization for dislocation
            IF (verbosity.GE.verbosity_max+1) THEN
                    WRITE(out,'(a)') 'Initialization for dislocation dipoles'
                    WRITE(out,*)
            END IF
            CALL RearrangeDDipole()
            CALL InitDDipole_ani(out)
    END IF

    ! ====Initialization for periodic calculation elastic calculation ======
    IF ( xImages.OR.yImages.OR.zImages ) THEN
            ! We need to correct displacement field, elastic strain, stress
            ! and energies to take into account periodic boundary conditions
            !   W. Cai, V. Bulatov, J. Chang, J. Li and S. Yip
            !   "Periodic image effects in dislocation modelling", Phil. Mag. 83 (2003), p. 539
            IF (.NOT.at_defined) THEN
                    WRITE(0,'(a)') 'Basis vectors at(1:3,1:3) undefined'
                    WRITE(0,'(a)') '  => could not apply periodic boundary conditions'
                    WRITE(0,'(a)') '     You need to define at(1:3,1:3) &
                        &or to set xImages, yImages and zImages to .false. '
                    STOP '< InitElasticity >'
            END IF
                    IF (verbosity.GE.verbosity_max+1) &
                            WRITE(out,'(a)') 'Calculate correction of conditional convergent sums &
                                &due to images of dislocations'
                    CALL Sum_Correction(epsilon_correction, sigma_correction, out)
    ELSE
            epsilon_correction = 0.d0
            sigma_correction = 0.d0
            IF (verbosity.GE.verbosity_max+1) &
                    WRITE(out,'(a)') 'No correction for periodicity added to &
                        &displacement and stress field'
    END IF

    ! ==== Cost function ===========================
    DO i=1, imDis
       uDis(:,i) = Array_Displacement(xpDis(:,i)) &
         + MatMul(eStrain+epsilon_correction, xpDis(:,i))
    END DO
    fCost = Sum( ( uDis(1:3,1:imDis) - u0Dis(1:3,1:imDis) )**2 )

    IF (verbosity.GE.verbosity_max+1) THEN
            WRITE(out,'(a,g20.6)') '  displacementCost = ', fCost
            WRITE(out,*)
    END IF

    IF (verbosity_backup.GT.0) verbosity = verbosity_backup

  END FUNCTION DisplacementCost

  SUBROUTINE Print_DisplacementCost(out)

    USE Disloc_elasticity_ani
    USE DDipoleModule
    IMPLICIT NONE

    INTEGER, intent(in) :: out

    WRITE(out,'(a)') ''
    IF (l_dislo) CALL PrintDislo(out)
    IF (l_DDipole) CALL PrintDDipole(out)
    WRITE(out,'(a)') ''

  END SUBROUTINE Print_DisplacementCost

END MODULE fitDisplacementModule
