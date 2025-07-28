MODULE EulerModule

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE EulerDefects(out)
    ! Apply Eulerian scheme to determine defect coordinates and periodicity
    ! vectors: x = x0 + u(x), where u(x) is the displacement calculated in x
    ! (eulerian coordinates)
    
    USE Babel_data
    USE PeriodRelaxModule
    USE PeriodModule
    USE PeachKoehlerModule
    USE Rearrange
    USE Disloc_elasticity_ani
    USE lineCouple_elasticity_ani
    USE LineForce_elasticity_ani
    USE DDipoleModule
    USE LoopModule
    IMPLICIT NONE
    INTEGER, intent(in) :: out

    REAL(kind(0.d0)), dimension(1:3,1:3) :: eStrain0, epsi, new_at
    INTEGER :: n, n_Euler
    LOGICAL :: converged

    IF (l_DDipole) THEN
            WRITE(0,'(a)') 'Eulerian coordinates not implemented with dislocation dipoles'
            STOP '< EulerDefects >'
    END IF
    IF (l_loop) THEN
            WRITE(0,'(a)') 'Eulerian coordinates not implemented with dislocation loops'
            STOP '< EulerDefects >'
    END IF

    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
            WRITE(out,'(a)') 'EULER: determine Eulerian coordinates of lattice defects and periodicity vectors'
            WRITE(out,'(a,i0,a)')  '    maximal number of iterations: max_Euler=', max_Euler, ')'
            WRITE(out,'(a,g14.6)') '    maximal tolerance for displacement convergency:&
                & delta_Euler = ', delta_Euler
            WRITE(out,*)
            DO n=1, nd
               WRITE(out,'(a,i0,a,3g14.6)') '  initial coordinates for dislo ', n,' : ', c0Dislo(1:3,n)
            END DO
            DO n=1, nlf
               WRITE(out,'(a,i0,a,3g14.6)') '  initial coordinates for line-force ', n,' : ', c0LineForce(1:3,n)
            END DO
            DO n=1, nlc
               WRITE(out,'(a,i0,a,3g14.6)') '  initial coordinates for line-force couple ', n,' : ', c0LineCouple(1:3,n)
            END DO
            WRITE(out,'(a)') '  initial unit cell vectors:'
            WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at0(1:3,1)
            WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at0(1:3,2)
            WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at0(1:3,3)
            WRITE(out,*)
    END IF
    ! Backup initial homogeneous strain
    eStrain0(:,:)=eStrain

    ! Initialization and checks
    IF ( ( xImages.OR.yImages.OR.zImages ) .AND. (.NOT.at_defined) ) THEN
            WRITE(0,'(a)') 'Basis vectors at(1:3,1:3) undefined'
            WRITE(0,'(a)') '  => could not apply periodic boundary conditions'
            WRITE(0,'(a)') '     You need to define at(1:3,1:3) &
                &or to set xImages, yImages and zImages to .false. '
            STOP '< EulerDefects >'
    END IF
    epsilon_correction = 0.d0
    sigma_correction = 0.d0


    ! Self-consistency loop
    Euler_scheme: DO n_Euler=1, max_Euler

            converged=.TRUE.

            ! Homogeneous strain induced by relaxation associated with lattice defects
            IF (induced_homogeneous_strain) THEN
                    CALL RelaxPeriod(epsi)
                    eStrain = eStrain0 + epsi
            ELSE
                    eStrain = eStrain0
            END IF

            ! Initialization for periodic calculation elastic calculation
            IF ( xImages.OR.yImages.OR.zImages ) &
                    CALL Sum_Correction(epsilon_correction, sigma_correction)

            ! Displacements of line-defects
            DO n=1, nd
               uDislo(1:3,n) = Displacement_on_dislo(n)
               IF (Any( Abs(cDislo(:,n)-c0Dislo(:,n)-uDislo(:,n)).GT.delta_Euler)) &
                        converged=.FALSE.
            END DO
            DO n=1, nlf
               uLineForce(1:3,n) = Displacement_on_LineForce(n)
               IF (Any( Abs(cLineForce(:,n)-c0LineForce(:,n)-uLineForce(:,n)).GT.delta_Euler)) &
                        converged=.FALSE.
            END DO
            DO n=1, nlc
               uLineCouple(1:3,n) = Displacement_on_LineCouple(n)
               IF (Any( Abs(cLineCouple(:,n)-c0LineCouple(:,n)-uLineCouple(:,n)).GT.delta_Euler)) &
                        converged=.FALSE.
            END DO

            ! New coordinates of line-defects
            cDislo(:,1:nd) = c0Dislo(:,1:nd) + uDislo(:,1:nd)
            cLineForce(:,1:nlf) = c0LineForce(:,1:nlf) + uLineForce(:,1:nlf)
            cLineCouple(:,1:nlc) = c0LineCouple(:,1:nlc) + uLineCouple(:,1:nlc)

            ! Recalculate dislocation cut lines and rotation matrix
            IF (verbosity.GE.verbosity_max+2) THEN
                    CALL RearrangeDislo(out)
            ELSE
                    CALL RearrangeDislo()
            END IF
            ! Reinitialize elastic quantities for dislocations
            CALL InitDisloc_ani()

            ! New periodicity vectors
            new_at(1:3,1:3) = MatMul( matId + eStrain, at0 )
            IF (Any( Abs(new_at(:,:)-at(:,:)).GT.delta_Euler)) converged=.FALSE.
            at(:,:) = new_at(:,:)

            IF (verbosity.GE.verbosity_max+1) THEN
                    WRITE(out,*)
                    WRITE(out,'(a,i0)') "Eulerian scheme: iter = ", n_Euler
                    WRITE(out,*)
                    DO n=1, nd
                       WRITE(out,'(a,i0,a,3g14.6)') '  new coordinates for dislo ', n,' : ', cDislo(1:3,n)
                    END DO
                    DO n=1, nlf
                       WRITE(out,'(a,i0,a,3g14.6)') '  new coordinates for line-force ', n,' : ', cLineForce(1:3,n)
                    END DO
                    DO n=1, nlc
                       WRITE(out,'(a,i0,a,3g14.6)') '  new coordinates for line-force couple ', n,' : ', cLineCouple(1:3,n)
                    END DO
                    WRITE(out,'(a)') '  new unit cell vectors:'
                    WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
                    WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
                    WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
                    WRITE(out,*)
                    WRITE(out,'(a)') '-------------'
                    WRITE(out,*)
            END IF

            IF (Converged) exit ! Self-consistency has been reached

    END DO Euler_scheme

    IF (verbosity.GE.verbosity_max) THEN
            IF (Converged) THEN
                    WRITE(out,'(a)') 'EULER: self-consistency has been reached'
            ELSE
                    WRITE(out,'(a)') 'EULER: maximal number of iterations has been reached'
            END IF
            WRITE(out,*)
            DO n=1, nd
               WRITE(out,'(a,i0,a,3g14.6)') '  new coordinates for dislo ', n,' : ', cDislo(1:3,n)
            END DO
            DO n=1, nlf
               WRITE(out,'(a,i0,a,3g14.6)') '  new coordinates for line-force ', n,' : ', cLineForce(1:3,n)
            END DO
            DO n=1, nlc
               WRITE(out,'(a,i0,a,3g14.6)') '  new coordinates for line-force couple ', n,' : ', cLineCouple(1:3,n)
            END DO
            WRITE(out,'(a)') '  new unit cell vectors:'
            WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
            WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
            WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
    END IF

  END SUBROUTINE EulerDefects

END MODULE EulerModule
