MODULE fitVitekModule

  IMPLICIT NONE
  SAVE

  ! Burgers and line vectors for the dislocation, and projection direction for
  ! displacement
  REAL(kind(0.d0)), dimension(1:3), private :: bVitek, lVitek, pVitek
  REAL(kind(0.d0)), private :: bNormVitek

  ! Number of data in Vitek map (arrows)
  INTEGER, private :: nVitek
  ! Differential displacement along Burgers vector direction
  REAL(kind(0.d0)), dimension(:), allocatable, private :: dUij0, dUij
  ! Atom indexes corresponding to each arrow
  INTEGER, dimension(:,:), allocatable, private :: indexVitek

  ! Number of atoms
  INTEGER, private :: imVitek
  ! Coordinates of each atoms
  REAL(kind(0.d0)), dimension(:,:), allocatable, private :: xpVitek
  ! Displacement of each atom
  REAL(kind(0.d0)), dimension(:,:), allocatable, private :: uVitek
  ! Logical telling if corresponding atom is used in Vitek map
  LOGICAL, dimension(:), allocatable, private :: incVitek

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE ReadVitek(inp, threshold)
    ! Read Vitek map in file connected to unit inp
    ! Only points for which the differential displacement is greater than
    ! threshold are kept

    IMPLICIT NONE
    INTEGER, intent(in) :: inp  ! Unit connected to Vitek file
    REAL(kind(0.d0)), intent(in) :: threshold

    CHARACTER(len=1) :: aTemp
    REAL(kind(0.d0)) :: dU
    INTEGER :: n, i, j, io
    REAL(kind(0.d0)), dimension(1:3) :: Xi, Xj

    ! Count number of uncommented lines => number of data in Vitek map
    !   and look for number of atoms
    imVitek=0
    nVitek=0
    DO
       CALL Comment(inp)
       READ(inp,*,iostat=io) Xi(1:3), Xj(1:3), dU, i, j
       IF (io.NE.0) exit
       IF (Abs(dU).LT.threshold) Cycle
       IF (i.GT.imVitek) imVitek=i
       IF (j.GT.imVitek) imVitek=j
       nVitek=nVitek+1
    END DO
    REWIND(inp)

    ! Read first line: comment containing Burgers vector
    READ(inp,*) aTemp, bVitek(1:3)
    ! Read second line: comment containing line direction
    READ(inp,*) aTemp, lVitek(1:3)

    ! Displacement projected along Burgers vector direction 
    bNormVitek = Sqrt( Sum( bVitek(:)**2 ) )
    pVitek(:) = bVitek(:) / bNormVitek

    ! Table allocation and initialization
    ALLOCATE(dUij0(1:nVitek))
    ALLOCATE(indexVitek(1:2,1:nVitek))
    ALLOCATE(xpVitek(1:3,1:imVitek)) 
    ALLOCATE(incVitek(1:imVitek))
    incVitek(:)=.FALSE.

    ! Read atom positions and differential displacement
    DO n=1, nVitek
       CALL Comment(inp)
       READ(inp,*) Xi(1:3), Xj(1:3), dU, i, j
       IF (Abs(dU).LT.threshold) Cycle
       dUij0(n) = dU
       indexVitek(1:2,n) = (/ i,j /)
       xpVitek(:,i) = Xi(:)
       xpVitek(:,j) = Xj(:)
       incVitek(i)=.TRUE.
       incVitek(j)=.TRUE.
    END DO

  END SUBROUTINE ReadVitek

  SUBROUTINE WriteVitek(dUijLocal,out)

    IMPLICIT NONE
    INTEGER, intent(in) :: out  ! Unit connected to Vitek file
    REAL(kind(0.d0)), dimension(:), intent(in) :: dUijLocal

    INTEGER :: n, i, j

    ! Read first line: comment containing Burgers vector
    WRITE(out,'(a,3g14.6,a)') '#', bVitek(1:3), ' (Burgers vector)'
    ! Read second line: comment containing line direction
    WRITE(out,'(a,3g14.6,a)') '#', lVitek(1:3), ' (line direction)'
    WRITE(out,'(a)') '# 1-3: 1st atom x and y coordinates'
    WRITE(out,'(a)') '# 4-6: 2nd atom x and y coordinates'
    WRITE(out,'(a)') '#   7: displacement difference uz2-uz1 along b'
    WRITE(out,'(a)') '# 8-9: atom indexes i and j'

    ! Write atom positions and differential displacement
    DO n=1, nVitek
       i = indexVitek(1,n)
       j = indexVitek(2,n)
       WRITE(out,'(7g14.6,2(1x,i0))') xpVitek(1:3,i), xpVitek(1:3,j), dUijLocal(n), i, j
    END DO

  END SUBROUTINE WriteVitek

  SUBROUTINE FitVitek(out)

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
    ALLOCATE(dUij(1:nVitek))
    ALLOCATE(uVitek(1:3,1:imVitek))

    !!! ==== DEBUG =================================
    !!fCost = VitekCost(x, nDim)
    !!OPEN(file='vitek.fit',unit=60,action='write')
    !!CALL WriteVitek(dUij, 60)
    !!CLOSE(60)
    !!STOP

    ! ==== First call =====
    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Minimization of Vitek map'
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
    CALL powell(VitekCost,x,xi,nDim,nDim,ftol,iter,fCost, print_VitekCost)

    ! ==== Get result ==========
    CALL Fit_ReadX(x, nDim)
    ! ... and apply constraints
    CALL ApplyConstraint()

    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'End of minimization'
            CALL Fit_PrintX(x, nDim, out)
            fCost = VitekCost(x, nDim)
            WRITE(out,'(a,g20.12)') '  fcost = ', fcost
            WRITE(out,*)
    END IF

    DeAllocate(xi)
    DeAllocate(x)

    ! ==== Write obtained Vitek map in a file =======
    OPEN(file='vitek.fit',unit=60,action='write')
    CALL WriteVitek(dUij,60)
    CLOSE(60)
    OPEN(file='vitek.residual',unit=61,action='write')
    CALL WriteVitek(dUij0-dUij,61)
    CLOSE(61)
    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Obtained Vitek map written in file vitek.fit'
            WRITE(out,'(a)') 'Difference between input and fitted Vitek map written in file vitek.residual'
            WRITE(out,*)
    END IF

    ! ==== The End =====
    DEALLOCATE(dUij)
    DEALLOCATE(uVitek)

  END SUBROUTINE FitVitek

  FUNCTION VitekCost(x, nDim) RESULT(fCost)

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

    INTEGER :: i, j, n
    INTEGER :: out, verbosity_backup
    REAL(kind(0.d0)) :: diff

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
            WRITE(out,'(a)') 'VitekCost CALL'
            WRITE(out,'(a,i0)') '  nDim = ', nDim
            DO n=1, nDim
               WRITE(out,'(a,i0,a,g20.12)') '  x(',n,') = ', x(n)
            END DO
            WRITE(out,*)
    END IF


    ! ==== Initialize dislocation definition (positions, Burgers)
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
    DO i=1, imVitek     ! Loop on atoms to calculate displacement
       IF (incVitek(i)) uVitek(:,i) = Array_Displacement(xpVitek(:,i)) &
         + MatMul(eStrain+epsilon_correction, xpVitek(:,i))
    END DO

    fCost = 0.d0
    DO n=1, nVitek
       i = indexVitek(1,n)
       j = indexVitek(2,n)
       

       dUij(n) = Dot_Product( uVitek(:,j) - uVitek(:,i), pVitek(:) )
       dUij(n) = dUij(n) - bNormVitek*aNint(dUij(n)/bNormVitek)
       diff = dUij(n) - dUij0(n)
       IF (Abs(diff+bNormVitek).LT.Abs(diff)) THEN
               diff = diff + bNormVitek
               dUij(n) = dUij(n) + bNormVitek
       ELSE IF (Abs(diff-bNormVitek).LT.Abs(diff)) THEN
               diff = diff - bNormVitek
               dUij(n) = dUij(n) - bNormVitek
       END IF
       fCost = fCost + diff**2
    END DO

    IF (verbosity.GE.verbosity_max+1) THEN
            WRITE(out,'(a,g20.6)') '  VitekCost = ', fCost
            WRITE(out,*)
    END IF

    IF (verbosity_backup.GT.0) verbosity = verbosity_backup

  END FUNCTION VitekCost

  SUBROUTINE Print_VitekCost(out)

    USE Disloc_elasticity_ani
    USE DDipoleModule
    IMPLICIT NONE

    INTEGER, intent(in) :: out

    WRITE(out,'(a)') ''
    IF (l_dislo) CALL PrintDislo(out)
    IF (l_DDipole) CALL PrintDDipole(out)
    WRITE(out,'(a)') ''

  END SUBROUTINE Print_VitekCost

END MODULE fitVitekModule
