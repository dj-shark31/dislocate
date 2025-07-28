MODULE disloc_elasticity_ani
  !------------------------------------------------------------------------------------
  ! Anisotropic elasticity for dislocations
  !----------------------------------------
  !   Ref.: J.P. Hirth and J. Lothe, Theory of Dislocations 
  !             (Wiley, New York, 2nd ed., 1982), pp.423-483
  !         J.P. Hirth and J. Lothe, 
  !             Anisotropic elastic solutions for line defects in high-symmetry cases
  !             J. Appl. Phys. 44 (1973), p. 1029
  !         J.D. Eshelby, W.T. Read and W. Shockley,
  !             Anisotropic elasticity with applications to dislocation theory
  !             Acta Met. 1 (1953), p. 251
  !
  ! First you need to call InitDisloc_ani
  ! This assumes that the following parameters have been defined in module
  ! Babel_data:
  !   lDislo(1:3): dislo line vector (normalized)
  !   bDislo(1:3): Burgers vector (A)
  !   cDislo(1:3): point belonging to the dislo line (A)
  !   rotDislo(1:3,1:3): matrix for changing ref. frame from cartesian to dislo
  !   CVoigt(1:6,1:6): elastic constants in Voigt notation (GPA)
  ! 
  ! Then, you obtain the displacement and the stress fields created by a dislo
  ! in a point M by calling
  !     CALL Dislo_Displacement_ani(R(1:3),u)
  !     CALL Dislo_Stress_ani(R(1:3),sigma)
  ! where R(1:3) are the cartesian coordinates of M (A)
  ! in output, you get u(1:3)=(/ux,uy,uz/), the displacement vector (A)
  !    and sigma(1:6), the stress tensor in Voigt notation (GPa)
  

  PUBLIC :: InitDisloc_ani, Dislo_Displacement_ani, Dislo_Stress_ani
  
  !========================================================
  ! Definition of the dislocations read in file 'input.dat'
  LOGICAL :: l_dislo      ! True if one wants to create a dislocation
  INTEGER :: nd                         ! Number of dislocations in the simulation box
  INTEGER, parameter :: max_nd=200        ! Maximal value for nd
  REAL(kind(0.d0)), dimension(1:3,1:max_nd) :: lDislo    ! Line vector
  REAL(kind(0.d0)), dimension(1:3,1:max_nd) :: bDislo    ! Burgers vector (in A)
  REAL(kind(0.d0)), dimension(1:3,1:max_nd) :: cutDislo, cut0Dislo  ! Direction defining discontinuity
  REAL(kind(0.d0)), dimension(1:3,1:max_nd) :: cDislo, c0Dislo      ! Point belonging to the dislo line (in A)
  REAL(kind(0.d0)), dimension(1:3,1:max_nd) :: uDislo    !   (displacement -> Eulerian coordinates)
  REAL(kind(0.d0)), dimension(1:3,1:3,1:max_nd) :: rotDislo, inv_rotDislo, rot0Dislo, inv_rot0Dislo ! Matrix for changing ref. frame

  !---------------------------------------------------------
  ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
  !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
  !   and the equation det[a(i;j)]=0 has 6 solutions for the variable p: root(1:6) 
  COMPLEX(kind(0.d0)), dimension(1:6,1:max_nd), save, private :: rootDislo
  ! Matrices used to calculate displacement and stress field created by a dislo
  COMPLEX(kind(0.d0)), dimension(1:3,1:6,1:max_nd), save, private :: displacement_matrix
  COMPLEX(kind(0.d0)), dimension(1:6,1:6,1:max_nd), save, private :: stress_matrix
  ! Angular dependence of the energy
  REAL(kind(0.d0)), dimension(1:max_nd), save, public :: dislo_angular_energy
  ! Prelogarithmic energy factor used to calculate dislocation elastic energy
  REAL(kind(0.d0)), dimension(1:max_nd), save, public :: dislo_prelog_energy_factor
  ! Matrix used to calculate interaction energy between dislocations
  COMPLEX(kind(0.d0)), dimension(1:3,1:6,1:max_nd), save, private :: interaction

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE InitDisloc_ani(out)
    ! Initialization of all quantities used for anisotropic elastic calculation     
    !   for a dislocation

    USE babel_data
    USE elasticity_ani
    USE elasticity_Stroh
    IMPLICIT NONE
    INTEGER, intent(in), optional :: out  ! output unit

    ! Elastic constants in initial and dislo axes
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: elastic_C, elastic_Cp

    ! Burgers vector and point force in dislo axes
    REAL(kind(0.d0)), dimension(1:3) :: burgers, force
    REAL(kind(0.d0)) :: K0, b2

    ! Outputs of subroutine Build_Elastic_A and inputs of Build_Elastic_D
    COMPLEX(kind(0.d0)), dimension(1:3,1:6) :: elastic_A, elastic_L
    COMPLEX(kind(0.d0)), dimension(1:3,1:3,1:3,1:6) :: elastic_B
    
    ! Prelogarithmic tensor appearing in elastic energy: E = 1/2 bi Kij bj ln(R/rc)
    REAL(kind(0.d0)), dimension(1:3,1:3) :: Kij

    REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
    INTEGER :: n, i, j
    REAL(kind(0.d0)) :: norm, rnd
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt_noisy
    REAL(kind(0.d0)), dimension(1:3) :: u1, u2

    ! Add some noise to elastic constants
    norm = CVoigt_noise*sqrt( Sum( CVoigt(1:6,1:6)**2 )/36. )
    IF (norm.GT.0.d0) THEN
            CVoigt_noisy = CVoigt(:,:) + CVoigt_random( norm, out )
            ! Check if elastic constant matrix is symmetric and definite positive
            CALL Check_CVoigt(CVoigt_noisy)
            IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
                    WRITE(out,*)
                    WRITE(out,'(a)') 'A noise has been added to elastic constants for creating dislocations'
                    WRITE(out,'(a,g13.6)') '  absolute amplitude of the noise: ', norm
                    CALL Print_CVoigt(CVoigt_noisy, out)
            END IF
    ELSE
            CVoigt_noisy=CVoigt
    END IF

    IF ( isotropic_CVoigt(CVoigt_noisy) ) THEN
            WRITE(0,*)
            WRITE(0,'(a)') 'Elastic constants are isotropic'
            WRITE(0,'(a)') 'You need to add some noise to elastic &
                &constants so as to break symmetry'
            WRITE(0,'(a)') 'Add to your input file "CVoigt_noise=1.d-4"'
            WRITE(0,'(a)') '  (do not use a noise lower than 1.d-6)'
            STOP '< InitDisloc_ani >'
    END IF

    ! Initialization
    elastic_C=Unpack_CVoigt(CVoigt_noisy)

    loop_dislo: DO n=1, nd

            ! Rotate elastic constants in the dislocation axes
            elastic_Cp = Rotate_Elastic_Constants(elastic_C, rotDislo(1:3,1:3,n) )
            ! Rotate Burgers vector in dislo axes: burgers(1) = edge component (if non 0)
            !                                      burgers(2) = 0
            !                                      burgers(3) = screw component
            burgers(1:3) = MatMul( rotDislo(1:3,1:3,n),bDislo(1:3,n) )
            
            ! No line-force for a dislocation
            force(1:3) = 0.d0

            IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
                    WRITE(out,'(a,i0)') ' Initialization for dislocation ', n
                    WRITE(out,'(a,3g14.6)') '    Burgers vector in rotated axes (A): ', burgers(1:3)
                    WRITE(out,'(a,3g14.6)') '    line direction in rotated axes : ', MatMul( rotDislo(1:3,1:3,n), lDislo(1:3,n) )
                    WRITE(out,'(a,3g14.6)') '    cut direction in rotated axes : ', MatMul( rotDislo(1:3,1:3,n), cutDislo(1:3,n) )
                    WRITE(out,'(a)') '  elastic constants in rotated axes'
                    CALL Print_Elastic_Constants(elastic_Cp,out)
            END IF


            IF (Present(out)) THEN
                    CALL InitStroh(elastic_Cp, rootDislo(:,n), elastic_A, elastic_B, elastic_L, out)
                    CALL Build_DStroh(burgers, force, rootDislo(:,n), elastic_A, elastic_B, elastic_L, &
                        displacement_matrix(:,:,n), stress_matrix(:,:,n), &
                        dislo_angular_energy(n), dislo_prelog_energy_factor(n), interaction(:,:,n), out)
            ELSE
                    CALL InitStroh(elastic_Cp, rootDislo(:,n), elastic_A, elastic_B, elastic_L)
                    CALL Build_DStroh(burgers, force, rootDislo(:,n), elastic_A, elastic_B, elastic_L, &
                        displacement_matrix(:,:,n), stress_matrix(:,:,n), &
                        dislo_angular_energy(n), dislo_prelog_energy_factor(n), interaction(:,:,n))
            END IF

            ! Displacement in cartesian coordinate
            displacement_matrix(:,:,n) = MatMul( inv_rotDislo(:,:,n), displacement_matrix(:,:,n) )

            IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
                    b2 = Sum(bDislo(1:3,n)**2)

                    ! prelogarithmic energy factor: E = K0*ln(R/r)
                    WRITE(out,'(a)') '  prelogarithmic energy factor: E = K0*ln(R/r)'
                    WRITE(out,'(a,g20.12,a)') '     K0 = ', dislo_prelog_energy_factor(n), &
                        '  -  units = elastic constant * distance^2 (GPa.A^2)'
                    WRITE(out,'(a,g20.12,a)') '        = ', factorE*dislo_prelog_energy_factor(n), &
                        '  -  units = energy / distance (eV.A^-1)'
                    WRITE(out,'(a,g20.12,a)') '        = ', dislo_prelog_energy_factor(n)/b2*4*pi, &
                        '*b^2/4pi  -  units = elastic constant (GPa)'

                    ! angular dependence of the energy
                    WRITE(out,'(a)') '  angular dependence of the energy'
                    WRITE(out,'(a,g20.12,a)') '    E0 = ', dislo_angular_energy(n), &
                        '  -  units = elastic constant * distance^2 (GPa.A^2)'
                    WRITE(out,'(a,g20.12,a)') '       = ', factorE*dislo_angular_energy(n), &
                        '  -  units = energy / distance (eV.A^-1)'
                    WRITE(out,'(a,g20.12,a)') '       = ', dislo_angular_energy(n)/b2*4*pi, &
                        '*b^2/4pi  -  units = elastic constant (GPa)'

                    ! Prelogarithmic tensor appearing in elastic energy: E = 1/2 bi Kij bj ln(R/rc)
                    DO i=1, 3 ; DO j=1, 3
                        Kij(i,j) = aImag( Sum( elastic_L(i,1:3)*elastic_L(j,1:3) ) &
                                - Sum( elastic_L(i,4:6)*elastic_L(j,4:6) ) )
                    END DO ; END DO
                    Kij(:,:) = Kij(:,:)/(2.d0*pi)
                    Kij(:,:) = MatMul( Transpose( rotDislo(:,:,n) ), MatMul( &
                        Kij(:,:), rotDislo(:,:,n) ) )
                    WRITE(out,'(a)') '  Stroh prelogarithmic tensor appearing in interaction energy: E = bi^1 Kij bj^2 ln(R/r)  (original axes)'
                    WRITE(out,'(a,3g20.12,a)') '           | ', Kij(1,1:3), ' | '
                    WRITE(out,'(a,3g20.12,a)') '     Kij = | ', Kij(2,1:3), ' | &
                        & units = elastic constant (GPa)'
                    WRITE(out,'(a,3g20.12,a)') '           | ', Kij(3,1:3), ' | '
                    K0 = factorE*0.5d0*Sum( bDislo(:,n) * MatMul( Kij(:,:), bDislo(:,n) ) )
                    WRITE(out,'(a,g20.12,a)') '    check for dislo self energy: K0 ?= 1/2 bi Kij bj = ', K0, &
                        '  -  units = energy / distance (eV.A^-1)'
                    WRITE(out,*)


                    ! Limited expansion of the displacement around cut plane
                    u1(1:3) = Dble( MatMul(displacement_matrix(1:3,1:3,n), rootDislo(1:3,n) ) )
                    u2(1:3) = -0.5d0*Dble( MatMul(displacement_matrix(1:3,1:3,n), rootDislo(1:3,n)**2 ) )
                    WRITE(out,'(a)') '  limited expansion of the displacement around the cut plane:'
                    WRITE(out,'(a)') '     u(x,y) ~ u(x,0) + u1*y/x + u2*(y/x)^2'
                    WRITE(out,'(a,3g20.12)') '    u1(1:3) = ', u1(1:3)
                    WRITE(out,'(a,3g20.12)') '    u2(1:3) = ', u2(1:3)
                    WRITE(out,*)
            END IF

    END DO loop_dislo
            
  END SUBROUTINE InitDisloc_ani

  FUNCTION Dislo_Displacement_ani(R, n) RESULT(u)
    ! Calculate total displacement in point R(1:3) arising from the dislocation n
    !   ref.: Hirth and Lothe, p.445, Eq.13.91

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the dislocation
    REAL(kind(0.d0)), dimension(1:3) :: u

    REAL(kind(0.d0)), dimension(2) :: ds
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    
    ! Location of the point in the reference frame of the dislo
    ds(:) = MatMul( rotDislo(1:2,:,n), R(:)-cDislo(:,n))

    IF (Sum(ds(:)**2).LT.1.d-12) THEN
            WRITE(0,'(a,3g14.6)') 'Point R(1:3) = ', R(1:3)
            WRITE(0,'(a,3g14.6)') 'Dislocation center cDislo(1:3) = ', cDislo(1:3,n)
            WRITE(0,'(2(a,g14.6))') ' Distance from dislocation line: x = ', ds(1), '  - y = ', ds(2)
            WRITE(0,'(a)') '  => cannot calculate displacement'
            STOP "< Dislo_Displacement_ani >"
    END IF

    ! Displacement in the reference frame of the dislo
    vector(1:3) = log( ds(1) + rootDislo(1:3,n)*ds(2) )
    u(:) = Dble( MatMul(displacement_matrix(:,1:3,n), vector(:)) )

  END FUNCTION Dislo_Displacement_ani

  FUNCTION Dislo_Stress_ani(R, n) RESULT(sigma)
    ! Calculate stress in point R(1:3) arising from the dislocation n
    !   ref.: Hirth and Lothe, p.445, Eq.13.92

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the dislocation
    ! Stress in Voigt notation
    REAL(kind(0.d0)), dimension(1:6) :: sigma

    REAL(kind(0.d0)), dimension(1:3) :: dR
    REAL(kind(0.d0)) :: x, y
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    REAL(kind(0.d0)), dimension(1:6) :: sigma1
    REAL(kind(0.d0)), dimension(1:3,1:3) :: stress, stress1
    
    ! Location of the point in the reference frame of the dislo
    dR(1:3) = R(1:3)-cDislo(1:3,n)
    x = Sum( rotDislo(1,1:3,n)*dR(1:3) )
    y = Sum( rotDislo(2,1:3,n)*dR(1:3) )

    ! Stress in the reference frame of the dislo
    vector(1:3) = 1.d0/(x+rootDislo(1:3,n)*y)
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )

    ! Stress in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotDislo(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotDislo(1:3,1:3,n) ) )
    sigma(1) = stress(1,1)
    sigma(2) = stress(2,2)
    sigma(3) = stress(3,3)
    sigma(4) = 0.5d0*( stress(2,3) + stress(3,2) )
    sigma(5) = 0.5d0*( stress(1,3) + stress(3,1) )
    sigma(6) = 0.5d0*( stress(1,2) + stress(2,1) )

  END FUNCTION Dislo_Stress_ani

  FUNCTION Dislo_gradStrain_ani(R, n) RESULT(gradk_Eij)
    ! Calculate strain gradient in point R(1:3) arising from the dislocation n

    USE babel_data, ONLY : inv_CVoigt
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the dislocation
    ! Strain gradient: gradk_Eij(i,j,k) = d[epsi(i,j)]/d[x(k)]
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3) :: gradk_Eij

    REAL(kind(0.d0)), dimension(1:3) :: dR
    REAL(kind(0.d0)) :: x, y
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    REAL(kind(0.d0)), dimension(1:6) :: sigma, sigma1, dx_epsi, dy_epsi
    REAL(kind(0.d0)), dimension(1:3,1:3) :: stress, stress1, dx_strain, dy_strain
    
    ! Location of the point in the reference frame of the dislo
    dR(1:3) = R(1:3)-cDislo(1:3,n)
    x = Sum( rotDislo(1,1:3,n)*dR(1:3) )
    y = Sum( rotDislo(2,1:3,n)*dR(1:3) )

    ! x derivative of the stress in the reference frame of the dislo
    vector(1:3) = -1.d0/(x+rootDislo(1:3,n)*y)**2
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )
    ! ... in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotDislo(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotDislo(1:3,1:3,n) ) )
    sigma(1) = stress(1,1)
    sigma(2) = stress(2,2)
    sigma(3) = stress(3,3)
    sigma(4) = 0.5d0*( stress(2,3) + stress(3,2) )
    sigma(5) = 0.5d0*( stress(1,3) + stress(3,1) )
    sigma(6) = 0.5d0*( stress(1,2) + stress(2,1) )
    ! ... corresponding strain derivative
    dx_epsi(:) = MatMul( inv_CVoigt, sigma )
    dx_strain(:,:) = Reshape( (/ &
                dx_epsi(1),0.5d0*dx_epsi(6),0.5d0*dx_epsi(5), &
                0.5d0*dx_epsi(6),dx_epsi(2),0.5d0*dx_epsi(4), &
                0.5d0*dx_epsi(5),0.5d0*dx_epsi(4),dx_epsi(3) /), (/3,3/) )

    ! y derivative of the stress in the reference frame of the dislo
    vector(1:3) = -rootDislo(1:3,n)/(x+rootDislo(1:3,n)*y)**2
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )
    ! ... in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotDislo(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotDislo(1:3,1:3,n) ) )
    sigma(1) = stress(1,1)
    sigma(2) = stress(2,2)
    sigma(3) = stress(3,3)
    sigma(4) = 0.5d0*( stress(2,3) + stress(3,2) )
    sigma(5) = 0.5d0*( stress(1,3) + stress(3,1) )
    sigma(6) = 0.5d0*( stress(1,2) + stress(2,1) )
    ! ... corresponding strain derivative
    dy_epsi(:) = MatMul( inv_CVoigt, sigma )
    dy_strain(:,:) = Reshape( (/ &
                dy_epsi(1),0.5d0*dy_epsi(6),0.5d0*dy_epsi(5), &
                0.5d0*dy_epsi(6),dy_epsi(2),0.5d0*dy_epsi(4), &
                0.5d0*dy_epsi(5),0.5d0*dy_epsi(4),dy_epsi(3) /), (/3,3/) )

    ! Strain gradient in cartesian coordinates
    gradk_Eij(:,:,1) = inv_rotDislo(1,1,n)*dx_strain(:,:) &
                + inv_rotDislo(1,2,n)*dy_strain(:,:)
    gradk_Eij(:,:,2) = inv_rotDislo(2,1,n)*dx_strain(:,:) &
                + inv_rotDislo(2,2,n)*dy_strain(:,:)
    gradk_Eij(:,:,3) = inv_rotDislo(3,1,n)*dx_strain(:,:) &
                + inv_rotDislo(3,2,n)*dy_strain(:,:)

  END FUNCTION Dislo_gradStrain_ani

  FUNCTION Dislo_Dislo_Interaction_ani(n,m, shift) RESULT(E)
    ! Calculate interaction elastic energy between dislocations n and m
    !   dislocations have to be parallel

    USE babel_data, ONLY : verbosity, verbosity_debug
    IMPLICIT NONE
    INTEGER, intent(in) :: n, m 
    REAL(kind(0.d0)), dimension(1:3), intent(in), optional :: shift
    !!$INTEGER, intent(in) :: out
    REAL(kind(0.d0)) :: E
    
    REAL(kind(0.d0)), dimension(1:3) :: dR, Kn, bm
    REAL(kind(0.d0)) :: x, y
    !!$REAL(kind(0.d0)) :: c, s
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    COMPLEX(kind(0.d0)) :: KnRoot
    COMPLEX(kind(0.d0)), dimension(1:6) :: EnRoot
    INTEGER :: nRoot

    ! Position of dislo m in the reference frame of dislo n
    dR(:) = cDislo(:,m) - cDislo(:,n)
    IF (Present(shift)) dR(:) = dR(:) + shift(:)
    x = Sum( rotDislo(1,:,n)*dR(:) )
    y = Sum( rotDislo(2,:,n)*dR(:) )

    ! Rotate Burgers vector of dislo m in the reference frame of dislo n
    bm(:) = MatMul( rotDislo(:,:,n),bDislo(:,m) )

    ! Interaction energy
    vector(1:3) = log( x + rootDislo(1:3,n)*y )
    Kn(:) = Dble(MatMul(interaction(:,1:3,n),vector(1:3)))
    E = -Sum( Kn(:)*bm(:) )

    !!$IF (verbosity.GE.verbosity_debug) THEN
            !!$WRITE(out,*)
            !!$WRITE(out,'(2(a,i0))') 'interaction energy between dislo ', &
                !!$n, ' and ', m
            !!$DO nRoot=1, 6
               !!$WRITE(out,*)
               !!$KnRoot = -0.5d0*Sum(bm(:)*interaction(:,nRoot,n))
               !!$WRITE(out,'(a,i1,2(a,g12.6))') '  a=', nRoot," :           b_i L(i,a) L(j,a) b'_j = ", &
                !!$Dble(KnRoot), ' + i * ', aImag(KnRoot)
               !!$WRITE(out,'(2(a,g12.6))') '                                     p_a = ', &
                !!$Dble(rootDislo(nRoot,n)), ' + i * ', aImag(rootDislo(nRoot,n))
               !!$WRITE(out,'(2(a,g12.6))') '                               x + p_a y = ', &
                !!$Dble(x+rootDislo(nRoot,n)*y), ' + i * ', aImag(x+rootDislo(nRoot,n)*y)
               !!$WRITE(out,'(2(a,g12.6))') '                         ln( x + p_a y ) = ', &
                !!$Dble(log(x+rootDislo(nRoot,n)*y)), ' + i * ', aImag(log(x+rootDislo(nRoot,n)*y))
               !!$EnRoot(nRoot) = KnRoot*log(x+rootDislo(nRoot,n)*y)
               !!$WRITE(out,'(2(a,g12.6))') "  b_i L(i,a) L(j,a) b'_j ln( x + p_a y ) = ", &
                !!$Dble(EnRoot(nRoot)), ' + i * ', aImag(EnRoot(nRoot))
            !!$END DO
            !!$WRITE(out,'(2(a,g12.6))') " Sum = ", Dble(Sum(EnRoot(:))), ' + i * ', aImag(Sum(EnRoot(:)))
            !!$WRITE(out,*)
    !!$END IF

  END FUNCTION Dislo_Dislo_Interaction_ani

  FUNCTION Dipole_Dislo_Interaction_ani(n,m, p, shift) RESULT(E)
    ! Calculate interaction elastic energy between dislocation dipole n-m
    !   and dislocation p
    !   dislocations have to be parallel

    IMPLICIT NONE
    INTEGER, intent(in) :: n, m, p
    REAL(kind(0.d0)), dimension(1:3), intent(in), optional :: shift
    REAL(kind(0.d0)) :: E
    
    REAL(kind(0.d0)), dimension(1:3) :: dR, Kn, bp
    REAL(kind(0.d0)) :: x, y, Ax, Ay
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector

    ! Position of dislo p in the reference frame of dislo n
    dR(:) = cDislo(:,p) - cDislo(:,n)
    IF (Present(shift)) dR(:) = dR(:) + shift(:)
    x = Sum( rotDislo(1,:,n)*dR(:) )
    y = Sum( rotDislo(2,:,n)*dR(:) )

    ! Position of dislo m in the reference frame of dislo n
    dR(:) = cDislo(:,m) - cDislo(:,n)
    Ax = Sum( rotDislo(1,:,n)*dR(:) )
    !!$Ay = Sum( rotDislo(2,:,n)*dR(:) )
    !!$IF (Ay.NE.0.d0) THEN
            !!$WRITE(0,*) 'n=', n, ' m=',m, ' p=', p, ' Ay=', Ay
            !!$STOP
    !!$END IF

    ! Rotate Burgers vector of dislo p in the reference frame of dislo n
    bp(:) = MatMul( rotDislo(:,:,n),bDislo(:,p) )

    ! Interaction energy
    !!$vector(1:3) = log(1.d0 - ( Ax + rootDislo(1:3,n)*Ay )/( x + rootDislo(1:3,n)*y ) )
    vector(1:3) = log(1.d0 -  Ax/( x + rootDislo(1:3,n)*y ) )
    Kn(:) = Dble(MatMul(interaction(:,1:3,n),vector(1:3)))
    E = Sum( Kn(:)*bp(:) )

  END FUNCTION Dipole_Dislo_Interaction_ani

  FUNCTION Dipole_Dipole_Interaction_ani(n,m, p,q, shift) RESULT(E)
    ! Calculate interaction elastic energy between dislocation dipoles n-m
    !   and p-q
    !   dislocations have to be parallel

    IMPLICIT NONE
    INTEGER, intent(in) :: n, m, p, q
    REAL(kind(0.d0)), dimension(1:3), intent(in), optional :: shift
    REAL(kind(0.d0)) :: E
    
    REAL(kind(0.d0)), dimension(1:3) :: dR, Kn, bm
    REAL(kind(0.d0)) :: xp, yp, xq, yq, Ax, Ay
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector

    ! Position of dislo p in the reference frame of dislo n
    dR(:) = cDislo(:,p) - cDislo(:,n)
    IF (Present(shift)) dR(:) = dR(:) + shift(:)
    xp = Sum( rotDislo(1,:,n)*dR(:) )
    yp = Sum( rotDislo(2,:,n)*dR(:) )

    ! Position of dislo q in the reference frame of dislo n
    dR(:) = cDislo(:,q) - cDislo(:,n)
    IF (Present(shift)) dR(:) = dR(:) + shift(:)
    xq = Sum( rotDislo(1,:,n)*dR(:) )
    yq = Sum( rotDislo(2,:,n)*dR(:) )

    ! Position of dislo m in the reference frame of dislo n
    dR(:) = cDislo(:,m) - cDislo(:,n)
    Ax = Sum( rotDislo(1,:,n)*dR(:) )
    Ay = Sum( rotDislo(2,:,n)*dR(:) )

    ! Rotate Burgers vector of dislo p in the reference frame of dislo n
    bm(:) = MatMul( rotDislo(:,:,n),bDislo(:,p) )

    ! Interaction energy
    vector(1:3) = log( ( 1.d0 - ( Ax + rootDislo(1:3,n)*Ay )/( xp + rootDislo(1:3,n)*yp ) ) &
        / ( 1.d0 - ( Ax + rootDislo(1:3,n)*Ay )/( xq + rootDislo(1:3,n)*yq ) ) )
    Kn(:) = Dble(MatMul(interaction(:,1:3,n),vector(1:3)))
    E = Sum( Kn(:)*bm(:) )

  END FUNCTION Dipole_Dipole_Interaction_ani

  SUBROUTINE InitDislo()

    IMPLICIT NONE

    l_dislo = .FALSE.
    nd = 0
    lDislo(:,:) = 0.d0
    bDislo(:,:) = 0.d0
    cDislo(:,:) = 0.d0
    cutDislo(:,:) = 0.d0
    rotDislo(:,:,:) = 0.d0
    inv_rotDislo(:,:,:) = 0.d0

  END SUBROUTINE InitDislo

  SUBROUTINE ReadDislo(inp)

    IMPLICIT NONE
    INTEGER, intent(in) :: inp

    NAMELIST /dislo/ lDislo, bDislo, cDislo, cutDislo, nd

    CALL InitDislo()

    READ(inp,nml=dislo)
    l_dislo=.true.

    ! Check number of dislo
    IF (nd.LE.0) THEN
            CALL InitDislo()
            RETURN
    ELSEIF (nd.GT.max_nd) THEN
            WRITE(0,'(a,i0)') 'Maximal number of dislocations (max_nd): ', max_nd
            WRITE(0,'(a,i0)') 'Current number of dislocations (nd): ', nd
            WRITE(0,'(a)') 'you need to change value in disloc_ani.f90 and to re-compile'
            STOP '< ReadDislo >'
    END IF

  END SUBROUTINE ReadDislo

  SUBROUTINE PrintDislo(out)
    ! Print dislocation definition on output unit out

    IMPLICIT NONE
    INTEGER, intent(in) :: out

    INTEGER :: n, i

    IF (.NOT.l_Dislo) THEN
            WRITE(out,'(a)') ' No dislocation defined'
            WRITE(out,*)
            RETURN
    END IF

    DO n=1, nd
            WRITE(out,'(a,i0)') '  Dislocation ', n
            WRITE(out,'(a,3g14.6)') '    line vector:        ', lDislo(1:3,n)
            WRITE(out,'(a,3g14.6)') '    Burgers vector (A): ', bDislo(1:3,n)
            WRITE(out,'(a,3g14.6)') '    center (A):         ', cDislo(1:3,n)
            WRITE(out,'(a,3g14.6)') '    cutting direction:  ', cutDislo(1:3,n)
            WRITE(out,'(a)') '    corresponding rotation matrix:'
            DO i=1, 3
               WRITE(out,'(6x,3f18.4)') rotDislo(i,1:3,n)
            END DO
    END DO 
    WRITE(out,*)

  END SUBROUTINE PrintDislo

  SUBROUTINE ScaleDislo(a)
    ! Multiply distances by a for dislo definition

    IMPLICIT NONE
    REAL(kind(0.d0)), intent(in) :: a

    ! Multiply dislo Burgers vector, center, and velocity by lattice parameter
    bDislo(:,1:nd) = a*bDislo(:,1:nd)
    cDislo(:,1:nd) = a*cDislo(:,1:nd)

  END SUBROUTINE ScaleDislo

  SUBROUTINE RotateDislo(rot, out, verbose)
    ! Rotate dislocations according to rotation matrix rot(1:3,1:3)
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_rot
    
    inv_rot = Transpose( rot )

    lDislo(1:3,1:nd) = MatMul( rot(1:3,1:3), lDislo(1:3,1:nd) )
    bDislo(1:3,1:nd) = MatMul( rot(1:3,1:3), bDislo(1:3,1:nd) )
    cutDislo(1:3,1:nd) = MatMul( rot(1:3,1:3), cutDislo(1:3,1:nd) )
    cDislo(1:3,1:nd) = MatMul( rot(1:3,1:3), cDislo(1:3,1:nd) )
    DO n=1, nd
            ! Beware: rotDislo and inv_rotDislo are not rank-2 tensors 
            rotDislo(:,:,n) = MatMul( rotDislo(:,:,n), inv_rot(:,:) )
            inv_rotDislo(:,:,n) = MatMul( rot(:,:), inv_rotDislo(:,:,n) )
    END DO

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for dislocations after rotation'
            CALL PrintDislo(out)
    END IF

  END SUBROUTINE RotateDislo

  SUBROUTINE TranslateDislo(u, out, verbose)
    ! Add displacement u(1:3) to dislocation defininition
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: u
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n

    DO n=1, nd
       cDislo(1:3,n) = cDislo(1:3,n) + u(1:3)
    END DO

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for dislocations after translation'
            CALL PrintDislo(out)
    END IF

  END SUBROUTINE TranslateDislo

END MODULE disloc_elasticity_ani
