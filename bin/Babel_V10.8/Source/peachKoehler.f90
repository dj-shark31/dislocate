MODULE PeachKoehlerModule

CONTAINS
  
  SUBROUTINE Print_fields_on_lineDefects(out)
    ! Print on output unit out the displacement and the stress seen by all
    ! dislocations
    
    USE Babel_data, ONLY : factorE, inv_CVoigt
    USE lineCouple_elasticity_ani
    USE disloc_elasticity_ani
    USE DDipoleModule
    IMPLICIT NONE
    INTEGER, intent(in) :: out

    INTEGER :: n1, i
    REAL(kind(0.d0)), dimension(1:6) :: sigma
    REAL(kind(0.d0)), dimension(1:3) :: F, u
    ! Strain gradient: gradk_Eij(i,j,k) = d[epsi(i,j)]/d[x(k)]
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3) :: gradk_Eij

    IF ( (nd.LE.0).AND.(nlc.LE.0).AND.(nDDipole.LE.0) ) RETURN

    WRITE(out,'(a)') '   forces given in [Energy]/[Distance]^2 -> eV/A²'
    WRITE(out,'(a,g14.6)') '   factor used to convert [Stress]*[Distance]^3 in [Energy]: &
                & factorE = ', factorE
    WRITE(out,'(a,g14.6,a)') '   ( 1 GPa * 1 A^3 = ', 1.d9*1.d-30/1.602176530000000045d-19,' eV )'

    IF (nd.GT.0) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Displacement, applied strain, stress, and force on all dislocations'
            dislo_loop: DO n1=1, nd
               u = Displacement_on_dislo(n1)
               sigma = Stress_on_dislo(n1)
               F = PK_force_dislo(sigma, n1)
               WRITE(out,*)
               WRITE(out,'(a,i0,a)') '  dislocation ', n1, ':'
               WRITE(out,'(a,3g14.6,a)') '    displacement: u(1:3) = ', u(:), ' !!!! EXPERIMENTAL !!!!'
               WRITE(out,'(a,6g14.6)') '    strain: e(1:6) = ', MatMul(inv_CVoigt(:,:), sigma(:) )
               WRITE(out,'(a,6g14.6)') '    stress: s(1:6) = ', sigma(1:6)
               WRITE(out,'(a,3g14.6)') '    Peach-Koehler force: F(1:3) = ', factorE*F(:)
            END DO dislo_loop
            WRITE(out,*)
    END IF
    IF (nDDipole.GT.0) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Displacement, applied strain, stress, and force on all dislocation dipoles'
            DDipole_loop: DO n1=1, nDDipole
               DO i=1, 2
                   u = Displacement_on_DDipole(n1, i)
                   sigma = Stress_on_DDipole(n1, i)
                   F = PK_force_DDipole(sigma, n1, i)
                   WRITE(out,*)
                   WRITE(out,'(2(a,i0),a)') '  dislocation ', i, ' of dipole ', n1, ':'
                   WRITE(out,'(a,3g14.6,a)')  '    displacement: u(1:3) = ', u(:), ' !!!! EXPERIMENTAL !!!!'
                   WRITE(out,'(a,6g14.6)')  '    strain: e(1:6) = ', MatMul(inv_CVoigt(:,:), sigma(:) )
                   WRITE(out,'(a,6g14.6)')  '    stress: s(1:6) = ', sigma(1:6)
                   WRITE(out,'(a,3g14.6)')  '    Peach-Koehler force: F(1:3) = ', factorE*F(:)
                END DO
            END DO DDipole_loop
            WRITE(out,*)
    END IF
    IF (nlc.GT.0) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Displacement, applied strain and stress on all line-force couples'
            lCouple_loop: DO n1=1, nlc
               u = Displacement_on_linecouple(n1)
               sigma = Stress_on_lineCouple(n1)
               gradk_Eij = gradStrain_on_lineCouple(n1)
               F = lineCouple_force(gradk_Eij, n1)
               WRITE(out,*)
               WRITE(out,'(a,i0,a)') '  line-force couple ', n1, ':'
               WRITE(out,'(a,3g14.6,a)') '    displacement: u(1:3) = ', u, ' !!!! EXPERIMENTAL !!!!'
               WRITE(out,'(a,6g14.6)') '    strain: e(1:6) = ', MatMul(inv_CVoigt(:,:), sigma(:) )
               WRITE(out,'(a,6g14.6)') '    stress: s(1:6) = ', sigma(1:6)
               WRITE(out,'(a,3g14.6)') '    core field force: F(1:3) = ', factorE*F(:)
            END DO lCouple_loop
            WRITE(out,*)
    END IF
    WRITE(out,'(a)') '==========================='
    WRITE(out,*)

  END SUBROUTINE Print_fields_on_lineDefects

  FUNCTION PK_Force_dislo(sVoigt, n1) RESULT(F)
    ! Calculate Peach-Koehler force on dislo n1 arising from stress sVoigt(1:6)

    USE disloc_elasticity_ani
    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:6), intent(in) :: sVoigt  ! Applied stress
    INTEGER, intent(in) :: n1   ! Dislo index
    REAL(kind(0.d0)), dimension(1:3) :: F       ! Peach-Koehler force

    REAL(kind(0.d0)), dimension(1:3,1:3) :: sigma       ! Applied stress
    REAL(kind(0.d0)), dimension(1:3) :: temp

    sigma(:,:) = Reshape( (/ sVoigt(1), sVoigt(6), sVoigt(5), &
                             sVoigt(6), sVoigt(2), sVoigt(4), &
                             sVoigt(5), sVoigt(4), sVoigt(3) /) , (/3,3/) )

    temp(:) = MatMul( sigma(:,:), bDislo(:,n1) )
    F(:) = CrossProduct( temp(:), lDislo(:,n1) )

  END FUNCTION PK_Force_dislo

  FUNCTION PK_Force_DDipole(sVoigt, n1, i) RESULT(F)
    ! Calculate Peach-Koehler force on dislo i from dipole n1 arising from stress sVoigt(1:6)

    USE DDipoleModule
    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:6), intent(in) :: sVoigt  ! Applied stress
    INTEGER, intent(in) :: n1, i   ! Dislo and dipole index
    REAL(kind(0.d0)), dimension(1:3) :: F       ! Peach-Koehler force

    REAL(kind(0.d0)), dimension(1:3,1:3) :: sigma       ! Applied stress
    REAL(kind(0.d0)), dimension(1:3) :: temp

    sigma(:,:) = Reshape( (/ sVoigt(1), sVoigt(6), sVoigt(5), &
                             sVoigt(6), sVoigt(2), sVoigt(4), &
                             sVoigt(5), sVoigt(4), sVoigt(3) /) , (/3,3/) )

    temp(:) = MatMul( sigma(:,:), bDDipole(:,n1) )
    F(:) = CrossProduct( temp(:), lDDipole(:,n1) )

    IF (i.EQ.2) F(:)=-F(:)

  END FUNCTION PK_Force_DDipole

  FUNCTION Stress_on_dislo(n1) RESULT(s)
    ! Calculate stress on dislocation n1 arising from all other stress sources

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE Rearrange
    USE PeriodModule

    IMPLICIT NONE
    INTEGER, intent(in) :: n1
    REAL(kind(0.d0)), dimension(1:6) :: s
   
    INTEGER :: n
    REAL(kind(0.d0)), dimension(1:3) :: R
    REAL(kind(0.d0)) :: dR2

    ! Initialization
    s(1:6) = 0.d0

    ! Point where the stress has to be calculated
    R(1:3) = cDislo(1:3,n1)

    ! Stress due to other dislocations in primitive unit-cell
    dislo_loop: DO n=1, nd
       IF (n.EQ.n1) Cycle
       dR2 = Sum( ( cDislo(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       s(1:6) = s(1:6) + Dislo_Stress_ani(R(1:3), n)
    END DO dislo_loop

    ! Stress due to dislocation dipoles in primitive unit-cell
    DDipole_loop: DO n=1, nDDipole
       s(1:6) = s(1:6) + DDipole_Stress_ani(R(1:3), n)
    END DO DDipole_loop

    ! Stress due to line-forces in primitive unit-cell
    force_loop: DO n=1, nlf
       dR2 = Sum( ( cLineForce(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       s(1:6) = s(1:6) + LineForce_Stress_ani(R(1:3), n)
    END DO force_loop

    ! Stress due to line-force couples in primitive unit-cell
    couple_loop: DO n=1, nlc
       IF (dislo_lineCouple(n1,n)) Cycle        ! Associated line-defects
       dR2 = Sum( ( cLineCouple(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       s(1:6) = s(1:6) + LineCouple_Stress_ani(R(1:3), n)
    END DO couple_loop

    ! Add stress due to image defects
    s(1:6) = s(1:6) + Array_stress(R, only_images=.true.)

    ! Add correction for image periodicity
    s(1:6) = s(1:6) + sigma_correction(1:6)

    ! Add stress corresponding to homogeneous strain
    s(1:6) = s(1:6) + sStrainVoigt(1:6)

  END FUNCTION Stress_on_dislo

  FUNCTION Stress_on_DDipole(n1, i) RESULT(s)
    ! Calculate stress on dislocation i from dipole n1 arising from all other stress sources

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE Rearrange
    USE PeriodModule

    IMPLICIT NONE
    INTEGER, intent(in) :: n1, i
    REAL(kind(0.d0)), dimension(1:6) :: s
   
    INTEGER :: n, j
    REAL(kind(0.d0)), dimension(1:3) :: R

    ! Initialization
    s(1:6) = 0.d0

    ! Point where the stress has to be calculated
    SELECT CASE(i)
    CASE(1)
            R(1:3) = c1DDipole(1:3,n1)     
            j=2
    CASE(2)
            R(1:3) = c2DDipole(1:3,n1)     
            j=1
    END SELECT

    ! Stress due to other dislocations in primitive unit-cell
    dislo_loop: DO n=1, nd
       s(1:6) = s(1:6) + Dislo_Stress_ani(R(1:3), n)
    END DO dislo_loop

    ! Stress due to dislocation dipoles in primitive unit-cell
    DDipole_loop: DO n=1, nDDipole
       IF (n.EQ.n1) THEN
               s(1:6) = s(1:6) + DDipole_Stress_ani(R(1:3), n, onlyDislo=j)
       ELSE
               s(1:6) = s(1:6) + DDipole_Stress_ani(R(1:3), n)
       END IF
    END DO DDipole_loop

    ! Stress due to line-forces in primitive unit-cell
    force_loop: DO n=1, nlf
       s(1:6) = s(1:6) + LineForce_Stress_ani(R(1:3), n)
    END DO force_loop

    ! Stress due to line-force couples in primitive unit-cell
    couple_loop: DO n=1, nlc
       s(1:6) = s(1:6) + LineCouple_Stress_ani(R(1:3), n)
    END DO couple_loop

    ! Add stress due to image defects
    s(1:6) = s(1:6) + Array_stress(R, only_images=.true.)

    ! Add correction for image periodicity
    s(1:6) = s(1:6) + sigma_correction(1:6)

    ! Add stress corresponding to homogeneous strain
    s(1:6) = s(1:6) + sStrainVoigt(1:6)

  END FUNCTION Stress_on_DDipole

  FUNCTION Stress_on_lineCouple(n1) RESULT(s)
    ! Calculate stress on line-force couple n1 arising from all other stress sources

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE Rearrange
    USE PeriodModule

    IMPLICIT NONE
    INTEGER, intent(in) :: n1
    REAL(kind(0.d0)), dimension(1:6) :: s
   
    INTEGER :: n
    REAL(kind(0.d0)), dimension(1:3) :: R
    REAL(kind(0.d0)) :: dR2

    ! Initialization
    s(1:6) = 0.d0

    ! Point where the stress has to be calculated
    R(1:3) = cLineCouple(1:3,n1)

    ! Stress due to other dislocations in primitive unit-cell
    dislo_loop: DO n=1, nd
       IF (dislo_lineCouple(n,n1)) Cycle        ! Associated line-defects
       dR2 = Sum( ( cDislo(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       s(1:6) = s(1:6) + Dislo_Stress_ani(R(1:3), n)
    END DO dislo_loop

    ! Stress due to dislocation dipoles in primitive unit-cell
    DDipole_loop: DO n=1, nDDipole
       s(1:6) = s(1:6) + DDipole_Stress_ani(R(1:3), n)
    END DO DDipole_loop

    ! Stress due to line-forces in primitive unit-cell
    force_loop: DO n=1, nlf
       dR2 = Sum( ( cLineForce(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       s(1:6) = s(1:6) + LineForce_Stress_ani(R(1:3), n)
    END DO force_loop

    ! Stress due to line-force couples in primitive unit-cell
    couple_loop: DO n=1, nlc
       IF (n.EQ.n1) Cycle
       dR2 = Sum( ( cLineCouple(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       s(1:6) = s(1:6) + LineCouple_Stress_ani(R(1:3), n)
    END DO couple_loop

    ! Add stress due to image defects
    s(1:6) = s(1:6) + Array_stress(R, only_images=.true.)

    ! Add correction for image periodicity
    s(1:6) = s(1:6) + sigma_correction(1:6)

    ! Add stress corresponding to homogeneous strain
    s(1:6) = s(1:6) + sStrainVoigt(1:6)

  END FUNCTION Stress_on_lineCouple

  FUNCTION gradStrain_on_lineCouple(n1) RESULT(gradk_Eij)
    ! Calculate strain gradient on line-force couple n1 arising from all other stress sources

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE Rearrange

    IMPLICIT NONE
    INTEGER, intent(in) :: n1
    ! Strain gradient: gradk_Eij(i,j,k) = d[epsi(i,j)]/d[x(k)]
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3) :: gradk_Eij
   
    INTEGER :: n
    REAL(kind(0.d0)), dimension(1:3) :: R
    REAL(kind(0.d0)) :: dR2

    ! Initialization
    gradk_Eij(:,:,:) = 0.d0

    ! Point where the strain gradient has to be calculated
    R(1:3) = cLineCouple(1:3,n1)

    ! Strain gradient due to other dislocations in primitive unit-cell
    dislo_loop: DO n=1, nd
       IF (dislo_lineCouple(n,n1)) Cycle        ! Associated line-defects
       dR2 = Sum( ( cDislo(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       gradk_Eij = gradk_Eij + Dislo_gradStrain_ani(R(1:3), n)
    END DO dislo_loop

    ! Stress due to dislocation dipoles in primitive unit-cell
    DDipole_loop: DO n=1, nDDipole
       gradk_Eij = gradk_Eij + DDipole_gradStrain_ani(R(1:3), n)
    END DO DDipole_loop

    ! Stress due to line-forces in primitive unit-cell
    force_loop: DO n=1, nlf
       dR2 = Sum( ( cLineForce(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       gradk_Eij = gradk_Eij + LineForce_gradStrain_ani(R(1:3), n)
    END DO force_loop

    ! Stress due to line-force couples in primitive unit-cell
    couple_loop: DO n=1, nlc
       IF (n.EQ.n1) Cycle
       dR2 = Sum( ( cLineCouple(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       gradk_Eij = gradk_Eij + LineCouple_gradStrain_ani(R(1:3), n)
    END DO couple_loop

  END FUNCTION gradStrain_on_lineCouple

  FUNCTION lineCouple_force(gradk_Eij, n1) RESULT(F)
    ! Force on line-force couple n1 arising from strain gradient
    USE lineCouple_elasticity_ani
    IMPLICIT NONE

    INTEGER, intent(in) :: n1
    ! Strain gradient: gradk_Eij(i,j,k) = d[epsi(i,j)]/d[x(k)]
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3), intent(in) :: gradk_Eij
    REAL(kind(0.d0)), dimension(1:3) :: F

    REAL(kind(0.d0)), dimension(1:3,1:3,1:3) :: gradk_Eij_rotated
    INTEGER :: k

    DO k=1, 3
       ! Apply rotation to components i and j of gradk_Eij(i,j,k)
       gradk_Eij_rotated(:,:,k) = MatMul( rotLineCouple(:,:,n1), &
                MatMul( gradk_Eij(:,:,k), inv_rotLineCouple(:,:,n1) ) )
       ! Corresponding force
       F(k) = gradk_Eij_rotated(1,1,k)*mxLineCouple(n1) &
                + gradk_Eij_rotated(2,2,k)*myLineCouple(n1)
    END DO

  END FUNCTION lineCouple_force

  FUNCTION Displacement_on_dislo(n1) RESULT(u)
    ! Calculate displacement on dislocation n1 arising from all other stress sources

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE Rearrange
    USE PeriodModule

    IMPLICIT NONE
    INTEGER, intent(in) :: n1
    REAL(kind(0.d0)), dimension(1:3) :: u
   
    INTEGER :: n, ids, n2
    REAL(kind(0.d0)), dimension(1:3) :: R, R2, s, noise
    REAL(kind(0.d0)) :: scalar, dR2

    ! Initialization
    u(1:3) = 0.d0

    ! Point where the stress has to be calculated
    R(1:3) = cDislo(1:3,n1)

    ! Small noise to ensure that the point is not exactly on a cut-line
    Call Random_Number(s(1:3))
    noise(1:3) = distance_zero*s(1:3)

    ! System to which the dislocation belongs
    ids = disloSystem(n1)

    ! Check if dislocations have been grouped in dipole
    IF ( (bClosedSystem(ids)).AND.(bDipoleSystem(ids)) ) THEN
            IF (modulo(n1,2).EQ.0) THEN
                    n2=n1-1
            ELSE
                    n2=n1+1
            END IF
    ELSE
            n2=n1
    END IF
    R2(1:3) = cDislo(1:3,n2)

    ! Displacement due to other dislocations in primitive unit-cell
    dislo_loop: DO n=1, nd
       IF ( (n.EQ.n1).OR.(n.EQ.n2) ) Cycle
       dR2 = Sum( ( cDislo(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) THEN
               ! Ghost dislocation: one cannot define its displacement
               u(:) = 0.d0
               Exit
       END IF
       u(1:3) = u(1:3) + 0.5d0*( Dislo_Displacement_ani(R(:)+noise(:), n) &
                + Dislo_Displacement_ani(R(:)-noise(:), n) )
    END DO dislo_loop

    ! Displacement due to dislocation dipoles
    DDipole_loop: DO n=1, nDDipole
       u(1:3) = u(1:3) + 0.5d0*( Dislo_Displacement_ani(R(1:3)+noise(:), n) &
               + Dislo_Displacement_ani(R(1:3)-noise(:), n) )
    END DO DDipole_loop

    ! Displacement due to line-forces in primitive unit-cell
    force_loop: DO n=1, nlf
       dR2 = Sum( ( cLineForce(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       u(1:3) = u(1:3) + 0.5d0*( LineForce_Displacement_ani(R(:)+noise(:), n) &
                + LineForce_Displacement_ani(R(:)-noise(:), n) )
    END DO force_loop

    ! Displacement due to line-force couples in primitive unit-cell
    couple_loop: DO n=1, nlc
       IF (dislo_lineCouple(n1,n)) Cycle        ! Associated line-defects
       dR2 = Sum( ( cLineCouple(:,n) - R(:) )**2 )
       IF ( dR2.LE.distance_zero2 ) Cycle
       u(1:3) = u(1:3) + 0.5d0*( LineCouple_Displacement_ani(R(:)+noise(:), n) &
                + LineCouple_Displacement_ani(R(:)-noise(:), n) )
    END DO couple_loop

    ! Add displacement due to image defects
    !   + correction for image periodicity
    !   + displacement corresponding to homogeneous strain
    u(1:3) = u(1:3) + 0.5d0*( Array_displacement(R(:)+noise(:), only_images=.true.) &
                    + Array_displacement(R(:)-noise(:), only_images=.true.) ) &
                    + MatMul( epsilon_correction(:,:) + eStrain(:,:), R(:) )

    ! Remove displacement parallel to line direction
    scalar = Sum(u(1:3)*lDislo(1:3,n1))
    u(1:3) = u(1:3) - scalar*lDislo(1:3,n1)

  END FUNCTION Displacement_on_dislo

  FUNCTION Displacement_on_DDipole(n1, i) RESULT(u)
    ! Calculate displacement on dislocation i (1 or 2) of dipole n1 arising from all other stress sources

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE Rearrange
    USE PeriodModule

    IMPLICIT NONE
    INTEGER, intent(in) :: n1, i
    REAL(kind(0.d0)), dimension(1:3) :: u
   
    INTEGER :: n
    REAL(kind(0.d0)), dimension(1:3) :: R, s, noise
    REAL(kind(0.d0)) :: scalar

    ! Initialization
    u(1:3) = 0.d0

    ! Point where the stress has to be calculated
    SELECT CASE(i)
    CASE(1)
            R(1:3) = c1DDipole(1:3,n1)     
    CASE(2)
            R(1:3) = c2DDipole(1:3,n1)     
    END SELECT

    ! Small noise to ensure that the point is not exactly on a cut-line
    Call Random_Number(s(1:3))
    noise(1:3) = distance_zero*s(1:3)

    ! Displacement due to other dislocations in primitive unit-cell
    dislo_loop: DO n=1, nd
       u(1:3) = u(1:3) + 0.5d0*( Dislo_Displacement_ani(R(:)+noise(:), n) &
                + Dislo_Displacement_ani(R(:)-noise(:), n) )
    END DO dislo_loop

    ! Displacement due to dislocation dipoles
    DDipole_loop: DO n=1, nDDipole
       IF (n.EQ.n1) Cycle
       u(1:3) = u(1:3) + 0.5d0*( DDipole_Displacement_ani(R(1:3)+noise(:), n) &
               + DDipole_Displacement_ani(R(1:3)-noise(:), n) )
    END DO DDipole_loop

    ! Displacement due to line-forces in primitive unit-cell
    force_loop: DO n=1, nlf
       u(1:3) = u(1:3) + 0.5d0*( LineForce_Displacement_ani(R(:)+noise(:), n) &
                + LineForce_Displacement_ani(R(:)-noise(:), n) )
    END DO force_loop

    ! Displacement due to line-force couples in primitive unit-cell
    couple_loop: DO n=1, nlc
       u(1:3) = u(1:3) + 0.5d0*( LineCouple_Displacement_ani(R(:)+noise(:), n) &
                + LineCouple_Displacement_ani(R(:)-noise(:), n) )
    END DO couple_loop

    ! Add displacement due to image defects
    !   + correction for image periodicity
    !   + displacement corresponding to homogeneous strain
    u(1:3) = u(1:3) + 0.5d0*( Array_displacement(R(:)+noise(:), only_images=.true.) &
                    + Array_displacement(R(:)-noise(:), only_images=.true.) ) &
                    + MatMul( epsilon_correction(:,:) + eStrain(:,:), R(:) )

    ! Remove displacement parallel to line direction
    scalar = Sum(u(1:3)*lDislo(1:3,n1))
    u(1:3) = u(1:3) - scalar*lDislo(1:3,n1)

  END FUNCTION Displacement_on_DDipole

  FUNCTION Displacement_on_lineForce(n1) RESULT(u)
    ! Calculate displacement on lineForce n1 arising from all other stress sources

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE Rearrange
    USE PeriodModule

    IMPLICIT NONE
    INTEGER, intent(in) :: n1
    REAL(kind(0.d0)), dimension(1:3) :: u
   
    INTEGER :: n
    REAL(kind(0.d0)), dimension(1:3) :: R, s, noise
    REAL(kind(0.d0)) :: scalar

    ! Initialization
    u(1:3) = 0.d0

    ! Point where the stress has to be calculated
    R(1:3) = cLineForce(1:3,n1)

    ! Small noise to ensure that the point is not exactly on a cut-line
    Call Random_Number(s(1:3))
    noise(:) = distance_zero*s(:)

    ! Displacement due to other dislocations in primitive unit-cell
    dislo_loop: DO n=1, nd
       u(:) = u(:) + 0.5d0*( Dislo_Displacement_ani(R(:)+noise(:), n) &
                + Dislo_Displacement_ani(R(:)-noise(:), n) )
    END DO dislo_loop

    ! Displacement due to dislocation dipoles
    DDipole_loop: DO n=1, nDDipole
       u(1:3) = u(1:3) + 0.5d0*( Dislo_Displacement_ani(R(1:3)+noise(:), n) &
               + Dislo_Displacement_ani(R(1:3)-noise(:), n) )
    END DO DDipole_loop

    ! Displacement due to line-forces in primitive unit-cell
    force_loop: DO n=1, nlf
       IF (n.EQ.n1) Cycle
       u(:) = u(:) + 0.5d0*( LineForce_Displacement_ani(R(:)+noise(:), n) &
                + LineForce_Displacement_ani(R(:)-noise(:), n) )
    END DO force_loop

    ! Displacement due to line-force couples in primitive unit-cell
    couple_loop: DO n=1, nlc
       u(:) = u(:) + 0.5d0*( LineCouple_Displacement_ani(R(:)+noise(:), n) &
                + LineCouple_Displacement_ani(R(:)-noise(:), n) )
    END DO couple_loop

    ! Add displacement due to image defects
    !   + correction for image periodicity
    !   + displacement corresponding to homogeneous strain
    u(:) = u(:) + 0.5d0*( Array_displacement(R(:)+noise(:), only_images=.true.) &
                    + Array_displacement(R(:)-noise(:), only_images=.true.) ) &
                    + MatMul( epsilon_correction(:,:) + eStrain(:,:), R(:) )

    ! Remove displacement parallel to line direction
    scalar = Sum(u(1:3)*lLineForce(1:3,n1))
    u(1:3) = u(1:3) - scalar*lLineForce(1:3,n1)

  END FUNCTION Displacement_on_lineForce

  FUNCTION Displacement_on_lineCouple(n1) RESULT(u)
    ! Calculate displacement on dislocation n1 arising from all other stress sources

    USE Babel_data
    USE disloc_elasticity_ani
    USE DDipoleModule
    USE LineForce_elasticity_ani
    USE LineCouple_elasticity_ani
    USE Rearrange
    USE PeriodModule

    IMPLICIT NONE
    INTEGER, intent(in) :: n1
    REAL(kind(0.d0)), dimension(1:3) :: u
   
    INTEGER :: n
    REAL(kind(0.d0)), dimension(1:3) :: R, s, noise
    REAL(kind(0.d0)) :: scalar

    ! Initialization
    u(1:3) = 0.d0

    ! Point where the stress has to be calculated
    R(1:3) = cLineCouple(1:3,n1)

    ! Small noise to ensure that the point is not exactly on a cut-line
    Call Random_Number(s(1:3))
    noise(1:3) = distance_zero*s(1:3)

    ! Displacement due to other dislocations in primitive unit-cell
    dislo_loop: DO n=1, nd
       IF (dislo_lineCouple(n,n1)) Cycle        ! Associated line-defects
       u(:) = u(:) + 0.5d0*( Dislo_Displacement_ani(R(:)+noise(:), n) &
                + Dislo_Displacement_ani(R(:)-noise(:), n) )
    END DO dislo_loop

    ! Displacement due to dislocation dipoles
    DDipole_loop: DO n=1, nDDipole
       u(1:3) = u(1:3) + 0.5d0*( Dislo_Displacement_ani(R(1:3)+noise(:), n) &
               + Dislo_Displacement_ani(R(1:3)-noise(:), n) )
    END DO DDipole_loop

    ! Displacement due to line-forces in primitive unit-cell
    force_loop: DO n=1, nlf
       u(:) = u(:) + 0.5d0*( LineForce_Displacement_ani(R(:)+noise(:), n) &
                + LineForce_Displacement_ani(R(:)-noise(:), n) )
    END DO force_loop

    ! Displacement due to line-force couples in primitive unit-cell
    couple_loop: DO n=1, nlc
       IF (n.EQ.n1) Cycle
       u(:) = u(:) + 0.5d0*( LineCouple_Displacement_ani(R(:)+noise(:), n) &
                + LineCouple_Displacement_ani(R(:)-noise(:), n) )
    END DO couple_loop

    ! Add displacement due to image defects
    !   + correction for image periodicity
    !   + displacement corresponding to homogeneous strain
    u(:) = u(:) + 0.5d0*( Array_displacement(R(:)+noise(:), only_images=.true.) &
                    + Array_displacement(R(:)-noise(:), only_images=.true.) ) &
                    + MatMul( epsilon_correction(:,:) + eStrain(:,:), R(:) )

    ! Remove displacement parallel to line direction
    scalar = Sum(u(1:3)*lLineCouple(1:3,n1))
    u(1:3) = u(1:3) - scalar*lLineCouple(1:3,n1)

  END FUNCTION Displacement_on_lineCouple

END MODULE PeachKoehlerModule
