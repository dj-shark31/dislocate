MODULE PeriodModule

  INTEGER, parameter, private :: verbosity_max=5

  ! Strain and stress to add to correct errors due to periodic images
  !   because of elastic energy calculation, we need to keep track of each
  !   correction separately
  REAL(kind(0.d0)), dimension(1:3,1:3), public, save :: epsilon_correction
  REAL(kind(0.d0)), dimension(1:3,1:3), private, save :: epsilon_correction_DDipole, &
        epsilon_correction_d, &
        epsilon_correction_lf, &
        epsilon_correction_lc, &
        epsilon_correction_loop
  REAL(kind(0.d0)), dimension(1:6), public, save :: sigma_correction, &
        sigma_correction_d, &
        sigma_correction_DDipole, &
        sigma_correction_lf, &
        sigma_correction_lc, &
        sigma_correction_loop

CONTAINS

  SUBROUTINE InitSumCorrection(out) 

    USE Babel_data
    USE loopModule
    USE disloc_elasticity_ani
    USE lineCouple_elasticity_ani
    USE LineForce_elasticity_ani
    USE DDipoleModule
    IMPLICIT NONE
    INTEGER, intent(in), optional :: out ! Output unit

    INTEGER :: nd_temp, nDDipole_temp, nlf_temp, nlc_temp, nLoop_temp

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
                    STOP '< InitSumCorrection >'
            END IF
            nd_temp=nd ; nDDipole_temp=nDDipole ; nlf_temp=nlf ; nlc_temp=nlc ; nLoop_temp=nLoop
            nd=0       ; nDDipole=0            ; nlf=0        ; nlc=0        ; nLoop=0
            ! Correction for loops
            IF (nLoop_temp.GT.0) THEN
                    nLoop=nLoop_temp
                    IF (verbosity.GE.verbosity_max) &
                            WRITE(out,'(a)') 'Calculate correction of conditional convergent sums &
                                &due to images of loops'
                    CALL Sum_Correction(epsilon_correction_loop, sigma_correction_loop, out)
                    nLoop=0
            END IF
            ! Correction for dislocations
            IF (nd_temp.GT.0) THEN
                    nd=nd_temp
                    IF (verbosity.GE.verbosity_max) &
                            WRITE(out,'(a)') 'Calculate correction of conditional convergent sums &
                                &due to images of dislocations'
                    CALL Sum_Correction(epsilon_correction_d, sigma_correction_d, out)
                    nd=0
            END IF
            ! Correction for dislocation dipoles
            IF (nDDipole_temp.GT.0) THEN
                    nDDipole=nDDipole_temp
                    IF (verbosity.GE.verbosity_max) &
                            WRITE(out,'(a)') 'Calculate correction of conditional convergent sums &
                                &due to images of dislocation dipoles'
                    CALL Sum_Correction(epsilon_correction_DDipole, sigma_correction_DDipole, out)
                    nDDipole=0
            END IF
            ! Correction for line-forces
            IF (nlf_temp.GT.0) THEN
                    nlf=nlf_temp
                    IF (verbosity.GE.verbosity_max) &
                            WRITE(out,'(a)') 'Calculate correction of conditional convergent sums &
                                &due to images of line-forces'
                    CALL Sum_Correction(epsilon_correction_lf, sigma_correction_lf, out)
                    nlf=0
            END IF
            ! Correction for line-force couples
            IF (nlc_temp.GT.0) THEN
                    nlc=nlc_temp
                    IF (verbosity.GE.verbosity_max) &
                            WRITE(out,'(a)') 'Calculate correction of conditional convergent sums &
                                &due to images of line-force couples'
                    CALL Sum_Correction(epsilon_correction_lc, sigma_correction_lc, out)
                    nlc=0
            END IF
            nd=nd_temp ; nDDipole=nDDipole_temp ; nlf=nlf_temp ; nlc=nlc_temp ; nLoop=nLoop_temp
            epsilon_correction = epsilon_correction_d &
                + epsilon_correction_DDipole &
                + epsilon_correction_lf &
                + epsilon_correction_lc &
                + epsilon_correction_loop
            sigma_correction = sigma_correction_d &
                + sigma_correction_DDipole &
                + sigma_correction_lf &
                + sigma_correction_lc &
                + sigma_correction_loop
    ELSE
            epsilon_correction = 0.d0
            sigma_correction = 0.d0
            IF (verbosity.GE.verbosity_max) &
                    WRITE(out,'(a)') 'No correction for periodicity added to &
                        &displacement and stress field'
    END IF

  END SUBROUTINE InitSumCorrection

  SUBROUTINE Sum_Correction(eErr, sErr, out)
    ! Correction of displacement field, elastic strain, stress
    ! and energies to take into account periodic boundary conditions
    !   W. Cai, V. Bulatov, J. Chang, J. Li and S. Yip
    !   "Periodic image effects in dislocation modelling", Phil. Mag. 83 (2003), p. 539

    USE Babel_data
    USE loopModule
    USE disloc_elasticity_ani
    USE LineCouple_elasticity_ani
    USE LineForce_elasticity_ani
    USE DDipoleModule
    USE Math
    IMPLICIT NONE
    ! Strain and stress to add to apply homogeneous strain 
    !   and to correct errors due to periodic images
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: eErr
    REAL(kind(0.d0)), dimension(1:6), intent(out) :: sErr
    INTEGER, intent(in), optional :: out ! Output unit

    REAL(kind(0.d0)), dimension(1:3) :: R0, s, u0, R
    REAL(kind(0.d0)), dimension(1:3,1:3) :: u, inv_at

    ! Point where calculation is performed
    IF (nLoop.GT.0) THEN
            r0(:) = Sum(cLoop(:,1:nLoop),2)/dble(nLoop) ! Loop gravity center
    ELSE IF (nd.GT.0) THEN
            R0(1:3) = Sum(cDislo(1:3,1:nd),2)/dble(nd)    ! Dislocation gravity center
    ELSE IF (nDDipole.GT.0) THEN
            R0(1:3) = ( Sum(c1DDipole(1:3,1:nDDipole),2)+Sum(c2DDipole(1:3,1:nDDipole),2) )/dble(2*nDDipole)
                    ! Dislocation dipole gravity center
    ELSE IF (nlf.GT.0) THEN
            R0(1:3) = Sum(cLineForce(1:3,1:nlf),2)/dble(nlf)    ! Line-force gravity center
    ELSE IF (nlc.GT.0) THEN
            R0(1:3) = Sum(cLineCouple(1:3,1:nlc),2)/dble(nlc)   ! Line-force gravity center
    ELSE
            R0(1:3) = 0.d0
    END IF
    Call Random_Number(s(1:3))
    R0(1:3) = R0(1:3) - 0.5d0*( s(1)*at(1:3,1) + s(2)*at(1:3,2) + s(3)*at(1:3,3) )
    ! Corresponding displacement
    u0(1:3) = Array_Displacement(R0) 

    ! Periodicity along x
    IF (xImages) THEN
            R(1:3) = R0(1:3) + at(1:3,1)
            u(1:3,1) = Array_Displacement(R) - u0(1:3)
    ELSE
            u(1:3,1) = 0.d0
    END IF

    ! Periodicity along y
    IF (yImages) THEN
            R(1:3) = R0(1:3) + at(1:3,2)
            u(1:3,2) = Array_Displacement(R) - u0(1:3)
    ELSE
            u(1:3,2) = 0.d0
    END IF

    ! Periodicity along z
    IF (zImages) THEN
            R(1:3) = R0(1:3) + at(1:3,3)
            u(1:3,3) = Array_Displacement(R) - u0(1:3)
    ELSE
            u(1:3,3) = 0.d0
    END IF

    ! Strain to add to correct errors due to periodic images
    CALL MatInv(at, inv_at)
    eErr(1:3,1:3) = -MatMul( u, inv_at )
    sErr(1:6) = MatMul( CVoigt, (/ eErr(1,1), &
                eErr(2,2), eErr(3,3), &
                eErr(2,3)+eErr(3,2), &
                eErr(1,3)+eErr(3,1), &
                eErr(1,2)+eErr(2,1) /) )

    IF (verbosity.LT.verbosity_max) RETURN
    IF (.NOT.Present(out)) RETURN
    WRITE(out,'(a,3g14.6)') '  Correction for periodicity calculated in R0 = ', R0(1:3)
    WRITE(out,'(a)') '  Displacement gradient to add for correcting periodicity'
    WRITE(out,'(a,3g14.6,a)') '     |', eErr(1,1:3), ' |'
    WRITE(out,'(a,3g14.6,a)') '     |', eErr(2,1:3), ' |'
    WRITE(out,'(a,3g14.6,a)') '     |', eErr(3,1:3), ' |'
    WRITE(out,'(a)') '  Stress tensor to add for correcting periodicity (GPa)'
    WRITE(out,'(a,3g14.6,a)') '     |', sErr(1), sErr(6), sErr(5), ' |'
    WRITE(out,'(a,3g14.6,a)') '     |', sErr(6), sErr(2), sErr(4), ' |'
    WRITE(out,'(a,3g14.6,a)') '     |', sErr(5), sErr(4), sErr(3), ' |'
    WRITE(out,'(a,g14.6,a)') '    corresponding pressure, P = ', -Sum(sErr(1:3))/3.d0, ' GPa'
    WRITE(out,'(a)') '  Displacement after correction:'
    WRITE(out,'(a,3g14.6)') '           in R0:             ', Array_Displacement(R0) + &
        MatMul( eErr, R0 )
    WRITE(out,'(a,3g14.6)') '           in R0 + at(1:3,1): ', Array_Displacement(R0+at(1:3,1)) + &
        MatMul( eErr, R0+at(:,1) )
    WRITE(out,'(a,3g14.6)') '           in R0 + at(1:3,2): ', Array_Displacement(R0+at(1:3,2)) + &
        MatMul( eErr, R0+at(:,2) )
    WRITE(out,'(a,3g14.6)') '           in R0 + at(1:3,3): ', Array_Displacement(R0+at(1:3,3)) + &
        MatMul( eErr, R0+at(:,3) )
    WRITE(out,'(a)') '  Stress after correction:'
    WRITE(out,'(a,6g14.6)') '           in R0:             ', Array_Stress(R0) + sErr(1:6)
    WRITE(out,'(a,6g14.6)') '           in R0 + at(1:3,1): ', Array_Stress(R0+at(1:3,1)) + sErr(1:6)
    WRITE(out,'(a,6g14.6)') '           in R0 + at(1:3,2): ', Array_Stress(R0+at(1:3,2)) + sErr(1:6)
    WRITE(out,'(a,6g14.6)') '           in R0 + at(1:3,3): ', Array_Stress(R0+at(1:3,3)) + sErr(1:6)
    WRITE(out,*)

  END SUBROUTINE Sum_Correction

  FUNCTION Array_Displacement(R0, only_images) RESULT(u)
    ! Calculate elastic displacement in point R(1:3) arising from all dislocations,
    !   all line-force couples and their images
    ! If only_images is present and set to true, 
    !   then only image defects are considered

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineCouple_elasticity_ani
    USE LineForce_elasticity_ani
    USE LoopModule
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R0        ! (A)
    LOGICAL, intent(in), optional :: only_images
    REAL(kind(0.d0)), dimension(1:3) :: u

    REAL(kind(0.d0)), dimension(1:3) :: R
    INTEGER :: n, ix, iy, iz

    ! Initialization
    u(1:3) = 0.d0
    ! Loop on all image cells
    images_loop: DO ix=-nxImages, nxImages
      DO iy=-nyImages, nyImages
         DO iz=-nzImages, nzImages
           IF (Present(only_images)) THEN
                   IF (only_images.AND.(ix.EQ.0).AND.(iy.EQ.0).AND.(iz.EQ.0))&
                        Cycle
           END IF
           R(1:3) = R0(1:3) + dble(ix)*at(1:3,1) + dble(iy)*at(1:3,2) &
                + dble(iz)*at(1:3,3)
           ! Loop on all dislocation loops
           DO n=1, nLoop
              u(:) = u(:) + Loops_displacement(R(:), n)
           END DO
           ! Loop on all dislocations
           DO n=1, nd
              u(:) = u(:) + Dislo_Displacement_ani(R(1:3), n)
           END DO
           ! Loop on all dislocation dipoles
           DO n=1, nDDipole
              u(:) = u(:) + DDipole_Displacement_ani(R(1:3), n)
           END DO
           ! Loop on all line-force
           DO n=1, nlf
              u(1:3) = u(1:3) + LineForce_Displacement_ani(R(1:3), n)
           END DO
           ! Loop on all line-force couples
           DO n=1, nlc
              u(1:3) = u(1:3) + LineCouple_Displacement_ani(R(1:3), n)
           END DO
    END DO ; END DO ; END DO images_loop

  END FUNCTION Array_Displacement

  FUNCTION Array_Stress(R0, only_images) RESULT(s)
    ! Calculate stress in point R(1:3) arising from all dislocations,
    !   line-force couples and their images
    ! If only_images is present and set to true, 
    !   then only image defects are considered

    USE babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R0        ! (A)
    LOGICAL, intent(in), optional :: only_images
    ! Stress in Voigt notation
    REAL(kind(0.d0)), dimension(1:6) :: s
    REAL(kind(0.d0)), dimension(1:3) :: R
    INTEGER :: n, ix, iy, iz

    s(1:6) = 0.d0
    ! Loop on all image cells
    images_loop: DO ix=-nxImages, nxImages
      DO iy=-nyImages, nyImages
         DO iz=-nzImages, nzImages
           IF (Present(only_images)) THEN
                   IF (only_images.AND.(ix.EQ.0).AND.(iy.EQ.0).AND.(iz.EQ.0))&
                        Cycle
           END IF
           R(1:3) = R0(1:3) + dble(ix)*at(1:3,1) + dble(iy)*at(1:3,2) &
                + dble(iz)*at(1:3,3)
           ! Loop on all dislocations to create dislocations
           DO n=1, nd
              s(1:6) = s(1:6) + Dislo_Stress_ani(R(1:3), n)        ! Stress
           END DO
           ! Loop on all dislocation dipoles to create dislocations
           DO n=1, nDDipole
              s(1:6) = s(1:6) + DDipole_Stress_ani(R(1:3), n)        ! Stress
           END DO
           ! Loop on all line-forces
           DO n=1, nlf
              s(1:6) = s(1:6) + LineForce_Stress_ani(R(1:3), n)        ! Stress
           END DO
           ! Loop on all line-force couples
           DO n=1, nlc
              s(1:6) = s(1:6) + LineCouple_Stress_ani(R(1:3), n)        ! Stress
           END DO
    END DO ; END DO ; END DO images_loop

  END FUNCTION Array_Stress

END MODULE PeriodModule
