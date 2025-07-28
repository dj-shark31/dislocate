MODULE lineForce_elasticity_ani
  !------------------------------------------------------------------------------------
  ! Anisotropic elasticity for line-force (this is the same as for dislocations)
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
  ! First you need to call InitLineForce_ani
  ! This assumes that the following parameters have been defined 
  !   lLineForce(1:3): lineForce line vector (normalized)
  !   fLineForce(1:3): force (in GPa.A)
  !   cLineForce(1:3): point belonging to the line (A)
  !   rotLineForce(1:3,1:3): matrix for changing ref. frame from cartesian 
  !   CVoigt(1:6,1:6): elastic constants in Voigt notation (GPA)
  ! 
  ! Then, you obtain the displacement and the stress fields created by a line-force
  ! in a point M by calling
  !     CALL LineForce_Displacement_ani(R(1:3),u)
  !     CALL LineForce_Stress_ani(R(1:3),sigma)
  ! where R(1:3) are the cartesian coordinates of M (A)
  ! in output, you get u(1:3)=(/ux,uy,uz/), the displacement vector (A)
  !    and sigma(1:6), the stress tensor in Voigt notation (GPa)
  

  PUBLIC :: InitLineForce_ani, LineForce_Displacement_ani, LineForce_Stress_ani
  
  !========================================================
  ! Definition of the line-forces
  LOGICAL :: LineForce      ! True if one wants to create a line-force
  INTEGER :: nlf         ! Number of line-forces
  INTEGER, parameter :: max_nlf=200        ! Maximal value for nlf
  REAL(kind(0.d0)), dimension(1:3,1:max_nlf) :: lLineForce  ! Line vector
  REAL(kind(0.d0)), dimension(1:3,1:max_nlf) :: fLineForce  ! Force (in GPa.A)
  REAL(kind(0.d0)), dimension(1:3,1:max_nlf) :: cutLineForce! Direction defining discontinuity
  REAL(kind(0.d0)), dimension(1:3,1:max_nlf) :: cLineForce  ! Point belonging to the line (in A)
  REAL(kind(0.d0)), dimension(1:3,1:max_nlf) :: c0LineForce !   (initial position -> Lagrangian coordinates)
  REAL(kind(0.d0)), dimension(1:3,1:max_nlf) :: uLineForce  !    (displacement -> Eulerian coordinates)
  REAL(kind(0.d0)), dimension(1:3,1:3,1:max_nlf) :: rotLineForce, inv_rotLineForce ! Matrix for changing ref. frame

  !---------------------------------------------------------
  ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
  !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
  !   and the equation det[a(i;j)]=0 has 6 solutions for the variable p: root(1:6) 
  COMPLEX(kind(0.d0)), dimension(1:6,1:max_nlf), save, private :: rootLineForce
  ! Matrices used to calculate displacement and stress field created by a line-force
  COMPLEX(kind(0.d0)), dimension(1:3,1:6,1:max_nlf), save, private :: displacement_matrix
  COMPLEX(kind(0.d0)), dimension(1:6,1:6,1:max_nlf), save, private :: stress_matrix
  ! Angular dependence of the energy
  REAL(kind(0.d0)), dimension(1:max_nlf), save, public :: lineForce_angular_energy
  ! Prelogarithmic energy factor used to calculate line-force elastic energy
  REAL(kind(0.d0)), dimension(1:max_nlf), save, public :: lineForce_prelog_energy_factor
  ! Matrix used to calculate interaction energy between dislocations
  COMPLEX(kind(0.d0)), dimension(1:3,1:6,1:max_nlf), save, private :: interaction

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE InitLineForce_ani(out)
    ! Initialization of all quantities used for anisotropic elastic calculation     
    !   for a line-force

    USE babel_data
    USE elasticity_ani
    USE elasticity_Stroh
    IMPLICIT NONE
    INTEGER, intent(in) :: out  ! output unit

    ! Elastic constants in initial and line-force axes
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: elastic_C, elastic_Cp

    ! Burgers vector and point force in line-force axes
    REAL(kind(0.d0)), dimension(1:3) :: burgers, force

    ! Outputs of subroutine Build_Elastic_A and inputs of Build_Elastic_D
    COMPLEX(kind(0.d0)), dimension(1:3,1:6) :: elastic_A, elastic_L
    COMPLEX(kind(0.d0)), dimension(1:3,1:3,1:3,1:6) :: elastic_B

    INTEGER :: i, j, n
    REAL(kind(0.d0)) :: norm, rnd
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt_noisy

    ! Add some noise to elastic constants
    norm = CVoigt_noise*sqrt( Sum( CVoigt(1:6,1:6)**2 )/36. )
    IF (norm.GT.0.d0) THEN
            DO i=1, 6
               DO j=1, i
                  CALL Random_Number(rnd)
                  CVoigt_noisy(i,j) = CVoigt(i,j) + norm*(0.5d0-rnd)
                  IF (i.NE.j) CVoigt_noisy(j,i) = CVoigt_noisy(i,j)
               END DO
            END DO
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,*)
                    WRITE(out,'(a)') 'A noise has been added to elastic constants for creating line-forces'
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
            STOP '< InitLineForce_ani >'
    END IF

    ! Initialization
    elastic_C=Unpack_CVoigt(CVoigt_noisy)

    loop_lineForce: DO n=1, nlf

            ! Rotate elastic constants in the line-force axes
            elastic_Cp = Rotate_Elastic_Constants(elastic_C, rotLineForce(1:3,1:3,n) )
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a,i0)') ' Initialization for line-force ', n
            END IF

            ! No Burgers vector for a line-force
            burgers(1:3) = 0.d0 
            
            ! Rotate line-force 
            force(1:3) = MatMul( rotLineForce(1:3,1:3,n),fLineForce(1:3,n) )

            CALL InitStroh(elastic_Cp, rootLineForce(1:6,n), elastic_A, elastic_B, elastic_L, out)
            CALL Build_DStroh(burgers, force, rootLineForce(:,n), elastic_A, elastic_B, elastic_L, &
                displacement_matrix(:,:,n), stress_matrix(:,:,n), &
                lineForce_angular_energy(n), lineForce_prelog_energy_factor(n), interaction(:,:,n), out)

            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a)') '  prelogarithmic energy factor: E = K*ln(R/r)'
                    WRITE(out,'(a,g20.12,a)') '     K = ', lineForce_prelog_energy_factor(n), &
                        '  -  units = elastic constant * distance^2 (GPa.A^2)'
                    WRITE(out,'(a)') '  angular dependence of the energy: '
                    WRITE(out,'(a,g20.12,a)') '    E0 = ', lineForce_angular_energy(n), &
                        '  -  units = elastic constant * distance^2 (GPa.A^2)'
                    WRITE(out,*)
            END IF

    END DO loop_lineForce
            
  END SUBROUTINE InitLineForce_ani

  FUNCTION LineForce_Displacement_ani(R, n) RESULT(u)
    ! Calculate total displacement in point R(1:3) arising from the line-force n
    !   ref.: Hirth and Lothe, p.445, Eq.13.91

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the line-force
    REAL(kind(0.d0)), dimension(1:3) :: u

    REAL(kind(0.d0)), dimension(1:3) :: dR, u1
    REAL(kind(0.d0)) :: x, y
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector, u_imag
    
    ! Location of the point in the reference frame of the line-force
    dR(1:3) = R(1:3)-cLineForce(1:3,n)
    x = Sum( rotLineForce(1,1:3,n)*dR(1:3) )
    y = Sum( rotLineForce(2,1:3,n)*dR(1:3) )

    IF (x*x+y*y.LT.1.d-12) THEN
            WRITE(0,'(a,3g14.6)') 'Point R(1:3) = ', R(1:3)
            WRITE(0,'(a,3g14.6)') 'Line-force center cLineForce(1:3) = ', cLineForce(1:3,n)
            WRITE(0,'(a,3g14.6)') '     dR(1:3) = ', dR(1:3)
            WRITE(0,'(2(a,g14.6))') ' Distance from line-force: x = ', x, '  - y = ', y
            WRITE(0,'(a)') '  => cannot calculate displacement'
            STOP "< LineForce_Displacement_ani >"
    END IF

    ! Displacement in the reference frame of the line-force
    vector(1:3) = log( x + rootLineForce(1:3,n)*y )
    u_imag(1:3) = MatMul(displacement_matrix(1:3,1:3,n), vector(1:3)) 
    u1(1:3) = Dble(u_imag(1:3))

    ! Displacement in cartesian coordinates
    u(1:3) = MatMul( inv_rotLineForce(1:3,1:3,n), u1(1:3) ) 

  END FUNCTION LineForce_Displacement_ani

  FUNCTION LineForce_Stress_ani(R, n) RESULT(sigma)
    ! Calculate stress in point R(1:3) arising from the line-force n
    !   ref.: Hirth and Lothe, p.445, Eq.13.92

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the line-force
    ! Stress in Voigt notation
    REAL(kind(0.d0)), dimension(1:6) :: sigma

    REAL(kind(0.d0)), dimension(1:3) :: dR
    REAL(kind(0.d0)) :: x, y
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    REAL(kind(0.d0)), dimension(1:6) :: sigma1
    REAL(kind(0.d0)), dimension(1:3,1:3) :: stress, stress1
    
    ! Location of the point in the reference frame of the line-force
    dR(1:3) = R(1:3)-cLineForce(1:3,n)
    x = Sum( rotLineForce(1,1:3,n)*dR(1:3) )
    y = Sum( rotLineForce(2,1:3,n)*dR(1:3) )

    ! Stress in the reference frame of the line-force
    vector(1:3) = 1.d0/(x+rootLineForce(1:3,n)*y)
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )

    ! Stress in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotLineForce(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotLineForce(1:3,1:3,n) ) )
    sigma(1) = stress(1,1)
    sigma(2) = stress(2,2)
    sigma(3) = stress(3,3)
    sigma(4) = 0.5d0*( stress(2,3) + stress(3,2) )
    sigma(5) = 0.5d0*( stress(1,3) + stress(3,1) )
    sigma(6) = 0.5d0*( stress(1,2) + stress(2,1) )

  END FUNCTION LineForce_Stress_ani

  FUNCTION LineForce_gradStrain_ani(R, n) RESULT(gradk_Eij)
    ! Calculate strain gradient in point R(1:3) arising from the line-force n

    USE babel_data, ONLY : inv_CVoigt
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the line-force
    ! Strain gradient: gradk_Eij(i,j,k) = d[epsi(i,j)]/d[x(k)]
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3) :: gradk_Eij

    REAL(kind(0.d0)), dimension(1:3) :: dR
    REAL(kind(0.d0)) :: x, y
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    REAL(kind(0.d0)), dimension(1:6) :: sigma, sigma1, dx_epsi, dy_epsi
    REAL(kind(0.d0)), dimension(1:3,1:3) :: stress, stress1, dx_strain, dy_strain
    
    ! Location of the point in the reference frame of the line-force
    dR(1:3) = R(1:3)-cLineForce(1:3,n)
    x = Sum( rotLineForce(1,1:3,n)*dR(1:3) )
    y = Sum( rotLineForce(2,1:3,n)*dR(1:3) )

    ! x derivative of the stress in the reference frame of the line-force
    vector(1:3) = -1.d0/(x+rootLineForce(1:3,n)*y)**2
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )
    ! ... in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotLineForce(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotLineForce(1:3,1:3,n) ) )
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

    ! y derivative of the stress in the reference frame of the line-force
    vector(1:3) = -rootLineForce(1:3,n)/(x+rootLineForce(1:3,n)*y)**2
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )
    ! ... in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotLineForce(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotLineForce(1:3,1:3,n) ) )
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
    gradk_Eij(:,:,1) = inv_rotLineForce(1,1,n)*dx_strain(:,:) &
                + inv_rotLineForce(1,2,n)*dy_strain(:,:)
    gradk_Eij(:,:,2) = inv_rotLineForce(2,1,n)*dx_strain(:,:) &
                + inv_rotLineForce(2,2,n)*dy_strain(:,:)
    gradk_Eij(:,:,3) = inv_rotLineForce(3,1,n)*dx_strain(:,:) &
                + inv_rotLineForce(3,2,n)*dy_strain(:,:)

  END FUNCTION LineForce_gradStrain_ani

  SUBROUTINE RotateLineForce(rot, out, verbose)
    ! Rotate line-forces according to rotation matrix rot(1:3,1:3)
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n, i
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_rot
    
    inv_rot = Transpose( rot )

    ! Rotate definition of the line-forces
    lLineForce(1:3,1:nlf) = MatMul( rot(1:3,1:3), lLineForce(1:3,1:nlf) )
    fLineForce(1:3,1:nlf) = MatMul( rot(1:3,1:3), fLineForce(1:3,1:nlf) )
    cutLineForce(1:3,1:nlf) = MatMul( rot(1:3,1:3), cutLineForce(1:3,1:nlf) )
    cLineForce(1:3,1:nlf) = MatMul( rot(1:3,1:3), cLineForce(1:3,1:nlf) )
    loop_force1: DO n=1, nlf
            ! Beware: rotLineForce and inv_rotLineForce are not rank-2 tensors 
            rotLineForce(1:3,1:3,n) = MatMul( rotLineForce(1:3,1:3,n), inv_rot(1:3,1:3) )
            inv_rotLineForce(1:3,1:3,n) = MatMul( rot(1:3,1:3), inv_rotLineForce(1:3,1:3,n) )
    END DO loop_force1

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for line-forces after rotation'
            CALL PrintLineForce(out)
    END IF

  END SUBROUTINE RotateLineForce

  SUBROUTINE TranslateLineForce(u, out, verbose)
    ! Add displacement u(1:3) to line-forces defininition
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: u
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n

    DO n=1, nlf
       cLineForce(1:3,n) = cLineForce(1:3,n) + u(1:3)
    END DO 

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for line-forces after translation'
            CALL PrintLineForce(out)
    END IF

  END SUBROUTINE TranslateLineForce

  SUBROUTINE PrintLineForce(out)
    ! Print dislocation dipole definition on output unit out

    IMPLICIT NONE
    INTEGER, intent(in) :: out

    INTEGER :: n, i

    IF (.NOT.LineForce) THEN
            WRITE(out,'(a)') ' No line-force defined'
            WRITE(out,*)
            RETURN
    END IF

    DO n=1, nlf
            WRITE(out,'(a,i0)') '  Line-force ', n
            WRITE(out,'(a,3g14.6)') '    line vector:        ', lLineForce(1:3,n)
            WRITE(out,'(a,3g14.6)') '    force (GPa.A):      ', fLineForce(1:3,n)
            WRITE(out,'(a,3g14.6)') '    center (A):         ', cLineForce(1:3,n)
            WRITE(out,'(a,3g14.6)') '    cutting direction:  ', cutLineForce(1:3,n)
            WRITE(out,'(a)') '    corresponding rotation matrix:'
            DO i=1, 3
               WRITE(out,'(6x,3f18.4)') rotLineForce(i,1:3,n)
            END DO
    END DO
    WRITE(out,*)

  END SUBROUTINE PrintLineForce

END MODULE LineForce_elasticity_ani
