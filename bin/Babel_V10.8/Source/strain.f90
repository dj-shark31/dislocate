MODULE strain_module

  IMPLICIT NONE
  ! stress(1:6,i): stress on atom
  ! pressure(i): pressure
  ! VonMises(i): equivalent Von-Misès shear stress
  REAL(kind(0.d0)), dimension(:,:), allocatable, save :: stress, elasticStrain
  REAL(kind(0.d0)), dimension(:), allocatable, save :: pressure, VonMises, Ebinding

  INTEGER, parameter, private :: verbosity_max=3

CONTAINS

  SUBROUTINE elastic_strain(uStrain, out)

    USE Babel_data
    USE disloc_elasticity_ani
    USE lineCouple_elasticity_ani
    USE LineForce_elasticity_ani
    USE DDipoleModule
    USE Math
    USE EnergyModule
    USE EnergyDipoleModule 
    USE PeriodModule
    USE PeachKoehlerModule
    USE Volterra_Module

    IMPLICIT NONE

    ! uStrain(1:3,i): total displacement of atom i
    REAL(kind(0.d0)), dimension(:,:), intent(out) :: uStrain
    INTEGER, intent(in) :: out

    INTEGER :: i, n, n_Euler
    REAL(kind(0.d0)), parameter :: UnTiers=1.d0/3.d0
    REAL(kind(0.d0)) :: Eel, cross
    ! Solid displacement that needs to be removed to keeped fixed the gravity center
    REAL(kind(0.d0)), dimension(1:3) :: uSolid, s, x, x0, u, dx
    REAL(kind(0.d0)), dimension(1:6) :: pVoigtImpurity
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
    REAL(kind(0.d0)) :: new_volume
    LOGICAL :: calculate_stress, calculate_displacement

    ! Check if stress needs to be calculated
    IF ( (out_stress).OR.(out_pressure).OR.(out_VonMises).OR.(out_Ebinding).OR.(out_elasticStrain) ) THEN
            calculate_stress=.TRUE.
    ELSE
            calculate_stress=.FALSE.
    END IF

    ! Check if displacement needs to be calculated
    IF (max_Euler.GE.1) THEN
            calculate_displacement=.TRUE.
    ELSE IF (outXyz.OR.outOnlyAtoms.OR.outCfg.OR.outGin.OR.outLisa.OR.outSiesta.OR.outNDM.OR.outLammps.OR.outPoscar) THEN
            IF (.NOT.initial) THEN
                    calculate_displacement=.TRUE.
            ELSE IF (out_displacement) THEN
                    calculate_displacement=.TRUE.
            ELSE
                    calculate_displacement=.FALSE.
            END IF
    ELSE
            calculate_displacement=.FALSE.
    END IF

    ! Allocation
    IF (Allocated(stress)) Deallocate(stress)
    IF (Allocated(pressure)) Deallocate(pressure)
    IF (Allocated(VonMises)) Deallocate(VonMises)
    IF (Allocated(elasticStrain)) Deallocate(elasticStrain)
    IF (Allocated(Ebinding)) Deallocate(Ebinding)
    uStrain(:,:)=0.d0
    IF (calculate_stress) THEN
            Allocate(stress(1:6,1:imm))
            stress(:,:)=0.d0
    END IF
    IF (out_pressure.OR.out_VonMises) THEN
            Allocate(pressure(1:imm))
            pressure(:)=0.d0
    END IF
    IF (out_VonMises) THEN
            Allocate(VonMises(1:imm))
            VonMises(:)=0.d0
    END IF
    IF (out_Ebinding) THEN
            ALLOCATE(Ebinding(1:imm))
            Ebinding(:)=0.d0
    END IF
    IF (out_elasticStrain) THEN
            ALLOCATE(elasticStrain(1:6,1:imm))
            elasticStrain(:,:)=0.d0
    END IF

    ! Initialization for elastic calculation
    CALL InitElasticity(out)

    ! Energy calculation
    IF ( (l_dislo.OR.l_lineCouple.OR.strain).AND..NOT.l_DDipole) THEN
            CALL ElasticEnergy(Eel, sigma_correction_d, sigma_correction_lc, out)
    ELSE IF (l_DDipole) THEN
            CALL ElasticEnergyDipole(Eel, sigma_correction_DDipole, out)
    END IF
    
    ! Calculate displacement and stress on each dislocation
    IF (verbosity.GE.verbosity_max) CALL Print_fields_on_lineDefects(out)

    IF (verbosity.GE.verbosity_max) THEN
            IF (l_dislo) WRITE(out,'(a)') 'Create dislocations'
            IF (l_DDipole) WRITE(out,'(a)') 'Create dislocation dipoles'
            IF (LineForce) WRITE(out,'(a)') 'Create line-forces'
            IF (l_LineCouple) WRITE(out,'(a)') 'Create line-force couples'
            IF (Strain) WRITE(out,'(a)') 'Apply homogeneous strain'
    END IF

    IF (calculate_displacement) THEN

            ! Loop on all atoms to create dislocations, line-forces and line-force couples
            DO i=1, im
               uStrain(:,i) = Array_Displacement(xp(:,i)) + MatMul(eStrain+epsilon_correction, xp(:,i))
            END DO
        
            ! Euler sef-consistency loop
            IF (max_Euler.GE.1) THEN
                    DO i=1, im
                       x0(:)=xp0(:,i) ; u(:)=uStrain(:,i) ; x(:)=x0(:)

                       DO n_Euler=1, max_Euler  ! Self-consistency loop for Eulerian coordinates
                          IF (ALL( Abs(x(:)-x0(:)-u(:)).LE.delta_Euler)) Exit
                          x(:) = x0(:) + u(:)
                          ! Check if the point crosses a dislo branch-cut and correct position
                          dx(:) = 0.d0
                          DO n=1, nd
                             cross = Branch_cutDislo_cross(x0(:), x(:), n)
                             dx(:) = dx(:) + cross*0.5d0*bDislo(:,n)
                          END DO
                          DO n=1, nDDipole
                             cross = Branch_cutDDipole_cross(x0(:), x(:), n)
                             dx(:) = dx(:) + cross*0.5d0*bDDipole(:,n)
                          END DO
                          ! Calculate displacement
                          u(:) = Array_Displacement(x(:)+dx(:)) + MatMul(eStrain+epsilon_correction, x(:)+dx(:))
                       END DO
                       IF ((n_Euler.GE.max_Euler).AND.(max_Euler.GT.0)) THEN
                               ! self-consistency could not been reached
                               IF (verbosity.GE.verbosity_max) THEN
                                       WRITE(out,'(a,i5,a,3g14.6)') 'EULER: atom ', i, ': xp0(1:3) = ', x0(1:3)
                                       WRITE(out,'(a,3g14.6)') '                    xp(1:3) = ', x(1:3)
                                       WRITE(out,'(a,3g14.6)') '                     u(1:3) = ', u(1:3)
                                       WRITE(out,'(a)') '  self-consistency could not been reached &
                                                &=> Lagrangian coordinates will be used for this atom'
                                       DO n=1, nd
                                          WRITE(out,'(a,i0,a,f0.0,a,3g14.6)') '  dislo ', n, &
                                                ' - cross = ', Branch_cutDislo_cross(x0(:), x(:), n), &
                                                ' - b(1:3) = ', bDislo(1:3,n)
                                       END DO
                                       DO n=1, nDDipole
                                          WRITE(out,'(a,i0,a,f0.0,a,3g14.6)') '  dislo dipole ', n, &
                                                ' - cross = ', Branch_cutDDipole_cross(x0(:), x(:), n), &
                                                ' - b(1:3) = ', bDDipole(1:3,n)
                                       END DO
                                       WRITE(out,*)
                               END IF
                               ! => Lagrangian coordinates
                               u(:) = Array_Displacement(x0(:)) + MatMul(eStrain+epsilon_correction, x0(:))
                               x(:) = x0(:) + u(:)
                       END IF
                       ! Final values
                       uStrain(:,i) = u(:)
                       xp(:,i) = x0(:) + u(:)
                    END DO
            END IF
    END IF

    ! Stress
    IF (calculate_stress) THEN
            IF (max_Euler.LT.1) THEN    ! Lagrangian coordinates
                    DO i=1, im
                       stress(1:6,i) = Array_Stress(xp0(1:3,i)) &           ! Stress
                                + sStrainVoigt(1:6) + sigma_correction(1:6)
                    END DO
            ELSE                        ! Eulerian coordinates
                    DO i=1, im
                       stress(1:6,i) = Array_Stress(xp(1:3,i)) &           ! Stress
                                + sStrainVoigt(1:6) + sigma_correction(1:6)
                    END DO
            END IF
    END IF

    ! Pressure
    IF (out_pressure.OR.out_VonMises) &
        pressure(1:im) = -UnTiers*( stress(1,1:im) + stress(2,1:im) + stress(3,1:im) ) 
    ! Von-Misès equivalent stress
    IF (out_VonMises) &
        VonMises(1:im) = Sqrt(1.5d0*( ( stress(1,1:im) + pressure(1:im) )**2 &
           + ( stress(2,1:im) + pressure(1:im) )**2 + ( stress(3,1:im) + pressure(1:im) )**2 &
           + 2.d0*(stress(4,1:im)**2 + stress(5,1:im)**2 + stress(6,1:im)**2 ) ))

    ! Elastic strain
    IF (out_elasticStrain) THEN
            elasticStrain(:,1:im) = MatMul(inv_CVoigt,stress(:,1:im))
    END IF

    ! Apply homogeneous strain to periodicity vectors
    IF (at_defined) THEN
            at(1:3,1:3) = MatMul( matId + eStrain, at0 )
            new_volume = abs(MatDet(at))
            IF (verbosity.GE.verbosity_max) THEN
                    IF (Strain) THEN
                            WRITE(out,'(a)') 'Homogeneous strain applied to unit cell vectors'
                            WRITE(out,'(a)') '  new unit cell vectors:'
                            WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
                            WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
                            WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
                            IF (im.LE.0) THEN
                                    WRITE(out,'(a,g20.12)')  '  corresponding volume ', new_volume
                            ELSE
                                    WRITE(out,'(2(a,g20.12),a)')  '  corresponding volume ', new_volume, &
                                        '  ( ', new_volume/dble(im), ' / atom )'
                            END IF
                            WRITE(out,*)
                    END IF
            END IF
    END IF

    ! Add displacement read in a file
    IF (read_uFile) THEN
            DO i=1, im
               uStrain(1:3,i) = uStrain(1:3,i) + uFile(1:3,i)
            END DO
    END IF

    ! Calculate solid displacement and remove it to keep fixed the gravity center
    IF (fixGravity .AND. (im.NE.0) ) THEN
            uSolid(1:3) = Sum( uStrain(1:3,1:im), 2)/dble(im)
            DO i=1, im
              uStrain(1:3,i) = uStrain(1:3,i) - uSolid(1:3)
            END DO
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a,3g14.6)') 'Solid displacement removed &
                        &to keep fixed the gravity center: u(1:3) = ', uSolid(1:3)
                    WRITE(out,*)
                    WRITE(out,'(a)') '==========================='
                    WRITE(out,*)
            END IF
    END IF

    !!$ IF (.NOT.initial) THEN
    !!$         ! initial=.FALSE. => apply displacement to atom coordinates
    !!$         xp(:,1:im) = xp0(:,1:im) + uStrain(:,1:im)              ! New cartesian coordinates
    !!$ ELSE
    !!$         xp(:,1:im) = xp0(:,1:im)
    !!$         at=at0
    !!$ END IF

    IF (out_Ebinding) THEN
            ! Calculate interaction energy between impurity and stress (Eshelby inclusion)
            pVoigtImpurity(1) = pImpurity(1,1)
            pVoigtImpurity(2) = pImpurity(2,2)
            pVoigtImpurity(3) = pImpurity(3,3)
            pVoigtImpurity(4) = 0.5d0*( pImpurity(2,3) + pImpurity(3,2) )
            pVoigtImpurity(5) = 0.5d0*( pImpurity(1,3) + pImpurity(3,1) )
            pVoigtImpurity(6) = 0.5d0*( pImpurity(1,2) + pImpurity(2,1) )
            DO i=1, im
               Ebinding(i) = Sum( pVoigtImpurity(:) * MatMul( inv_CVoigt(:,:), stress(:,i) ) )
            END DO
    END IF


  END SUBROUTINE elastic_strain

  SUBROUTINE InitElasticity(out)

    USE Babel_data
    USE Math
    USE Disloc_elasticity_ani
    USE DDipoleModule
    !USE LoopModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE PeriodModule
    USE PeriodRelaxModule
    USE EulerModule
    IMPLICIT NONE
    INTEGER, intent(in), optional :: out
    INTEGER :: n! nd_temp, nDDipole_temp, nlf_temp, nlc_temp, nLoop_temp
    REAL(kind(0.d0)), dimension(1:3,1:3) :: epsi, eStrain0

    INTEGER :: i

    IF (l_dislo) THEN             ! Initialization for dislocation
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a)') 'Initialization for dislocations'
                    WRITE(out,*)
            END IF
            CALL InitDisloc_ani(out)
    END IF
    IF (l_DDipole) THEN             ! Initialization for dislocation dipoles
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a)') 'Initialization for dislocation dipoles'
                    WRITE(out,*)
            END IF
            CALL InitDDipole_ani(out)
    END IF
    IF (lineForce) THEN         ! Initialization for line-forces
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a)') 'Initialization for line-forces'
                    WRITE(out,*)
            END IF
            CALL InitLineForce_ani(out)
    END IF
    IF (l_LineCouple) THEN        ! Initialization for line-force couples
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a)') 'Initialization for line-force couples'
                    WRITE(out,*)
            END IF
            CALL InitLineCouple_ani(out)
    END IF

    ! Save initial (Lagrangian) coordinates
    at0(:,:) = at(:,:)
    IF (l_dislo) THEN
            DO n=1, nd
               c0Dislo(:,n) = cDislo(:,n)
               cut0Dislo(:,n) = cutDislo(:,n)
               rot0Dislo(:,:,n) = rotDislo(:,:,n)
               inv_rot0Dislo(:,:,n) = inv_rotDislo(:,:,n)
            END DO
    END IF
    IF (lineForce) THEN
            DO n=1, nlf
               c0LineForce(:,n) = cLineForce(:,n)
            END DO
    END IF
    IF (l_lineCouple) THEN
            DO n=1, nlc
               c0LineCouple(:,n) = cLineCouple(:,n)
            END DO
    END IF

    ! Homogeneous strain
    IF (.NOT.strain) eStrain = 0.d0

    ! Backup initial homogeneous strain
    eStrain0(:,:)=eStrain
    
    !! Calculate Eulerian coordinates of lattice defects and periodicity vectors
    !IF (max_Euler.GT.0) CALL EulerDefects(out)

    ! Add homogeneous strain induced by line defects
    IF (induced_homogeneous_strain) THEN
            CALL RelaxPeriod(epsi,out)
            eStrain = eStrain0 + epsi
    ELSE
            eStrain = eStrain0
    END IF

    IF (matNorm2(eStrain).GT.1.d-6) strain=.true.
    ! Corresponding stress
    sStrainVoigt(1:6) = MatMul( CVoigt, (/ eStrain(1,1), eStrain(2,2), eStrain(3,3), &
        eStrain(2,3)+eStrain(3,2), eStrain(1,3)+eStrain(3,1), eStrain(1,2)+eStrain(2,1) /) )
    IF ((xImages.OR.yImages.OR.zImages).AND.(verbosity.GE.verbosity_max)) THEN
            WRITE(out,*)
            WRITE(out,'(a)') '--------------'
            WRITE(out,*)
            WRITE(out,'(a)') 'Homogeneous strain, applied and resulting&
                & from all line-defects'
            DO i=1,3
               WRITE(out,'(a,i1,a,3(g13.6,2x),a)') '  e(',i,',1:3) = | ', &
                    eStrain(i,1:3), ' |'
            END DO
            WRITE(out,*)
            WRITE(out,'(a)') '--------------'
            WRITE(out,*)
    END IF

    ! Initialization for periodic calculation elastic calculation
    ! We need to correct displacement field, elastic strain, stress
    ! and energies to take into account periodic boundary conditions
    !   W. Cai, V. Bulatov, J. Chang, J. Li and S. Yip
    !   "Periodic image effects in dislocation modelling", Phil. Mag. 83 (2003), p. 539
    CALL InitSumCorrection(out)

    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
    END IF

  END SUBROUTINE InitElasticity

END MODULE strain_module
