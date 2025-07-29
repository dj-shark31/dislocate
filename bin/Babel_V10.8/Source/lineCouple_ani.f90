MODULE LineCouple_elasticity_ani
  !------------------------------------------------------------------------------------
  ! Anisotropic elasticity for line force couples
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
  ! First you need to call InitCouple_ani
  ! This assumes that the following parameters have been defined in module
  !   lLineCouple(1:3): line vector (normalized)
  !   mxLineCouple, myLineCouple: force moments in x and y directions (GPa.A^2)
  !   xLineCouple, yLineCouple: corresponding x and y directions 
  !   cLineCouple(1:3): point belonging to the couple line (A)
  !   rotLineCouple(1:3,1:3): matrix for changing ref. frame from cartesian 
  !   CVoigt(1:6,1:6): elastic constants in Voigt notation (GPA)
  ! 
  ! Then, you obtain the displacement and the stress fields created 
  ! by a line-force couple in a point M by calling
  !     CALL LineCouple_Displacement_ani(R(1:3),u)
  !     CALL LineCouple_Stress_ani(R(1:3),sigma)
  ! where R(1:3) are the cartesian coordinates of M (A)
  ! in output, you get u(1:3)=(/ux,uy,uz/), the displacement vector (A)
  !    and sigma(1:6), the stress tensor in Voigt notation (GPa)

  PUBLIC :: InitLineCouple_ani, LineCouple_Displacement_ani, LineCouple_Stress_ani
  
  !========================================================
  ! Definition of the line-force couples read in file 'input.dat'
  LOGICAL :: l_lineCouple     ! True if one wants to create a line-force couple
  INTEGER :: nlc                         ! Number of line-force couples
  INTEGER, parameter :: max_nlc=200        ! Maximal value for nlc
  REAL(kind(0.d0)), dimension(1:3,1:max_nlc) :: lLineCouple  ! Line vector
  REAL(kind(0.d0)), dimension(1:max_nlc) :: mxLineCouple, myLineCouple  ! Amplitude of the force moments (in GPa.A^2)
  REAL(kind(0.d0)), dimension(1:3,1:max_nlc) :: xLineCouple, yLineCouple  ! Direction of the force moments (in GPa.A^2)
  REAL(kind(0.d0)), dimension(1:3,1:max_nlc) :: cLineCouple  ! Point belonging to the couple line (in A)
  REAL(kind(0.d0)), dimension(1:3,1:max_nlc) :: c0LineCouple !   (initial position -> Lagrangian coordinates)
  REAL(kind(0.d0)), dimension(1:3,1:max_nlc) :: uLineCouple  !    (displacement -> Eulerian coordinates)
  REAL(kind(0.d0)), dimension(1:3,1:3,1:max_nlc) :: rotLineCouple, inv_rotLineCouple ! Matrix for changing ref. frame

  !---------------------------------------------------------
  ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
  !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
  !   and the equation det[a(i;j)]=0 has 6 solutions for the variable p: root(1:6) 
  COMPLEX(kind(0.d0)), dimension(1:6,1:max_nlc), save, private :: rootLineCouple
  ! Matrices used to calculate displacement and stress field created
  COMPLEX(kind(0.d0)), dimension(1:3,1:6,1:max_nlc), save, private :: displacement_matrix
  COMPLEX(kind(0.d0)), dimension(1:6,1:6,1:max_nlc), save, private :: stress_matrix
  ! Energy factor used to calculate line-couple elastic energy
  REAL(kind(0.d0)), dimension(1:max_nlc), save, public :: lineCouple_energy_factor

  ! Matrix giving interaction energy with an external applied strain
  REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3,1:max_nlc), save, private :: inter_lineCouple
  REAL(kind(0.d0)), dimension(1:9,1:9,1:max_nlc), save :: inter_lineCouple_Voigt

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE InitLineCouple

    IMPLICIT NONE

    l_lineCouple=.FALSE.
    nlc=0
    lLineCouple(1:3,:) = 0.d0
    mxLineCouple(:) = 0.d0
    xLineCouple(1:3,:) = 0.d0
    myLineCouple(:) = 0.d0
    yLineCouple(1:3,:) = 0.d0
    cLineCouple(1:3,:) = 0.d0
    
  END SUBROUTINE InitLineCouple

  SUBROUTINE ReadLineCouple(inp)

    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    INTEGER :: nLineCouple 

    NAMELIST /lineCouple/ nLineCouple, lLineCouple, mxLineCouple, xLineCouple, myLineCouple, yLineCouple, cLineCouple

    CALL InitLineCouple()
    nLineCouple=0

    READ(inp,nml=lineCouple)
    nlc=nLineCouple
    l_lineCouple=.true.

    IF (nlc.LE.0) THEN
            WRITE(0,'(a)') 'You need to enter a positive number of line-force couples (nlc)'
            STOP '< ReadLineCouple >'
    ELSEIF (nlc.GT.max_nlc) THEN
            WRITE(0,'(a,i0)') 'Maximal number of line-force couples (max_nlc): ', max_nlc
            WRITE(0,'(a,i0)') 'Current number of line-force couples (nlc): ', nlc
            WRITE(0,'(a)') 'you need to change value in data.f90 and to re-compile'
            STOP '< ReadLineCouple >'
    END IF

  END SUBROUTINE ReadLineCouple

  SUBROUTINE PrintLineCouple(out)
    ! Print dislocation dipole definition on output unit out

    IMPLICIT NONE
    INTEGER, intent(in) :: out

    INTEGER :: n, i

    IF (.NOT.l_LineCouple) THEN
            WRITE(out,'(a)') ' No line-force dipole defined'
            WRITE(out,*)
            RETURN
    END IF

    loop_couple2: DO n=1, nlc
            WRITE(out,'(a,i0)') '  Line-force dipole ', n
            WRITE(out,'(a,3g14.6)') '    line vector:        ', lLineCouple(1:3,n)
            WRITE(out,'(a,g14.6,a,3g14.6)') '    force-moments (GPa.A^2): Mx = ', mxLineCouple(n), &
                ' in direction ', xLineCouple(1:3,n)
            WRITE(out,'(a,g14.6,a,3g14.6)') '                             My = ', myLineCouple(n), &
                ' in direction ', yLineCouple(1:3,n)
            WRITE(out,'(a,3g14.6)') '    center (A):         ', cLineCouple(1:3,n)
            WRITE(out,'(a)') '    corresponding rotation matrix:'
            DO i=1, 3
               WRITE(out,'(6x,3f18.4)') rotLineCouple(i,1:3,n)
            END DO
    END DO loop_couple2
    WRITE(out,*)

  END SUBROUTINE PrintLineCouple

  SUBROUTINE ScaleLineCouple(a)
    ! Multiply distances by a for line-force couples

    IMPLICIT NONE
    REAL(kind(0.d0)), intent(in) :: a

    ! Multiply center of line force couple by lattice parameter
    cLineCouple(:,1:nlc) = a*cLineCouple(:,1:nlc)

  END SUBROUTINE ScaleLineCouple

  SUBROUTINE InitLineCouple_ani(out)
    ! Initialization of all quantities used for anisotropic elastic calculation     
    !   for a line-force couple

    USE babel_data
    USE elasticity_ani
    USE elasticity_Stroh
    USE Math
    IMPLICIT NONE
    INTEGER, intent(in) :: out  ! output unit

    ! Elastic constants in initial and line-force axes
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: elastic_C, elastic_Cp

    ! Outputs of subroutine Build_Elastic_A and inputs of Build_Elastic_D
    COMPLEX(kind(0.d0)), dimension(1:3,1:6) :: elastic_A, elastic_L
    COMPLEX(kind(0.d0)), dimension(1:6) :: elastic_D
    COMPLEX(kind(0.d0)), dimension(1:3,1:3,1:3,1:6) :: elastic_B
    COMPLEX(kind(0.d0)) :: pn, pm, factor, Iy
    COMPLEX(kind(0.d0)), dimension(1:3) :: u
    COMPLEX(kind(0.d0)), dimension(1:6) :: sigma
    COMPLEX(kind(0.d0)), dimension(1:6,1:6) :: elastic_AL

    REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
    REAL(kind(0.d0)) :: nSign, mSign

    INTEGER :: i, j, k, l, m, n, p, nRoot, mRoot
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
                    WRITE(out,'(a)') 'A noise has been added to elastic constants for creating line-force couples'
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
            STOP '< InitLineCouple_ani >'
    END IF

    ! Initialization
    elastic_C=Unpack_CVoigt(CVoigt_noisy)

    loop_couple: DO n=1, nlc

            ! Rotate elastic constants in the line-force axes
            elastic_Cp = Rotate_Elastic_Constants(elastic_C, rotLineCouple(1:3,1:3,n) )
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a,i0)') ' Initialization for line-force couple ', n
            END IF

            CALL InitStroh(elastic_Cp, rootLineCouple(1:6,n), elastic_A, elastic_B, elastic_L, out)

            ! Matrix used to calculate displacement and stress fields
            factor=-1.d0/cmplx(0.d0,2.d0*pi,kind(0.d0))
            loop_roots: DO nRoot=1, 6
               pn = rootLineCouple(nRoot,n)
               IF (nRoot.EQ.4) factor=-factor
               elastic_D(nRoot) = -2.d0*( elastic_A(1,nRoot)*mxLineCouple(n) &
                        + elastic_A(2,nRoot)*pn*myLineCouple(n) )
               u(1:3) = factor*elastic_A(1:3,nRoot)*elastic_D(nRoot)
               sigma(1) = Sum( elastic_B(1,1,1:3,nRoot)*u(1:3) )
               sigma(2) = Sum( elastic_B(2,2,1:3,nRoot)*u(1:3) )
               sigma(3) = Sum( elastic_B(3,3,1:3,nRoot)*u(1:3) )
               sigma(4) = 0.5d0*Sum( (elastic_B(2,3,1:3,nRoot)+elastic_B(3,2,1:3,nRoot))*u(1:3) )
               sigma(5) = 0.5d0*Sum( (elastic_B(1,3,1:3,nRoot)+elastic_B(3,1,1:3,nRoot))*u(1:3) )
               sigma(6) = 0.5d0*Sum( (elastic_B(1,2,1:3,nRoot)+elastic_B(2,1,1:3,nRoot))*u(1:3) )
               displacement_matrix(1:3,nRoot,n) = - u(1:3)
               stress_matrix(1:6,nRoot,n) = sigma(1:6)
            END DO loop_roots

            ! Energy prefactor
            elastic_AL(:,:) = 0.d0
            factor = 0.d0
            nSign=+1.d0
            DO nRoot=1, 6
               IF (nRoot.EQ.4) nSign=-1.d0
               pn = rootLineCouple(nRoot,n)
               mSign=1.d0
               DO mRoot=1, 6
                  IF (mRoot.EQ.4) mSign=-1.d0
                  pm = rootLineCouple(mRoot,n)
                  IF ((aImag(pn).GE.0.d0).AND.(aImag(pm).LE.0.d0)) THEN
                          Iy = -(1.d0 + pn*pm)/(pn-pm)**2
                  ELSE IF ((aImag(pn).LE.0.d0).AND.(aImag(pm).GE.0.d0)) THEN
                          Iy = (1.d0 + pn*pm)/(pn-pm)**2
                  ELSE
                          Cycle
                  END IF
                  elastic_AL(nRoot,mRoot) = nSign*mSign*Iy &
                        * Sum( elastic_A(1:3,nRoot) * elastic_L(1:3,mRoot) )
               END DO
            END DO

            factor = Cmplx( 0.d0, 1.d0/(8.d0*pi) )*Sum(elastic_D(:) * &
                        MatMul(elastic_AL(:,:), elastic_D(:) ) )

            lineCouple_energy_factor(n) = Dble(factor)


            !--------------------------------------------------------
            ! Matrix giving interaction energy with an external strain
            inter_lineCouple(:,:,:,:,n) = &
                Build_Inter_lineCouple(elastic_Cp, rootLineCouple(1:6,n), elastic_A)

            ! ... in Voigt notation
            DO m=1,9
               CALL inv_Voigt_index(m,i,j)
               DO p=1,9
                  CALL inv_Voigt_index(p,k,l)
                  inter_lineCouple_Voigt(m,p, n) = inter_lineCouple(i,j,k,l, n)
               END DO
            END DO

            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a)') '  energy prefactor: E = -K/r^2'
                    WRITE(out,'(a,g20.12,a)') '     K = ', lineCouple_energy_factor(n), &
                        '  -  units = elastic constant * distance^4 (GPa.A^4)'
                    WRITE(out,'(a,g20.12,a)') '       = ', factorE*lineCouple_energy_factor(n), &
                        '  -  units = energy * distance (eV.A)'
                    WRITE(out,*)
                    WRITE(out,'(a)') '  interaction matrix with external strain'
                    WRITE(out,'(a)') '    (E = - e_ij L_ijmn M_mn for e_ij in line-couple axis)'
                    DO m=1,9
                       WRITE(out,'(a,9g14.6,a)') '   | ', inter_lineCouple_Voigt(m,1:9,n), ' |'
                    END DO
            END IF

    END DO loop_couple
            
  END SUBROUTINE InitLineCouple_ani

  FUNCTION LineCouple_Displacement_ani(R, n) RESULT(u)
    ! Calculate total displacement in point R(1:3) arising from line-force couple n

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the line-force couple
    REAL(kind(0.d0)), dimension(1:3) :: u

    REAL(kind(0.d0)), dimension(1:3) :: dR, u1
    REAL(kind(0.d0)) :: x, y
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector, u_imag
    
    ! Location of the point in the reference frame of the line-force couple
    dR(1:3) = R(1:3)-cLineCouple(1:3,n)
    x = Sum( rotLineCouple(1,1:3,n)*dR(1:3) )
    y = Sum( rotLineCouple(2,1:3,n)*dR(1:3) )

    IF (x*x+y*y.LT.1.d-12) THEN
            WRITE(0,'(a,3g14.6)') 'Point R(1:3) = ', R(1:3)
            WRITE(0,'(a,3g14.6)') 'Line-force couple center cCouple(1:3) = ', cLineCouple(1:3,n)
            WRITE(0,'(a,3g14.6)') '     dR(1:3) = ', dR(1:3)
            WRITE(0,'(2(a,g14.6))') ' Distance from line: x = ', x, '  - y = ', y
            WRITE(0,'(a)') '  => cannot calculate displacement'
            STOP "< LineCouple_Displacement_ani >"
    END IF

    ! Displacement in the reference frame of the line-force
    vector(1:3) = 1.d0/( x + rootLineCouple(1:3,n)*y )
    u_imag(1:3) = MatMul(displacement_matrix(1:3,1:3,n), vector(1:3)) 
    u1(1:3) = Dble(u_imag(1:3))

    ! Displacement in cartesian coordinates
    u(1:3) = MatMul( inv_rotLineCouple(1:3,1:3,n), u1(1:3) ) 

  END FUNCTION LineCouple_Displacement_ani

  FUNCTION LineCouple_Stress_ani(R, n) RESULT(sigma)
    ! Calculate stress in point R(1:3) arising from the line-force couple n

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the line-force couple
    ! Stress in Voigt notation
    REAL(kind(0.d0)), dimension(1:6) :: sigma

    REAL(kind(0.d0)), dimension(1:3) :: dR
    REAL(kind(0.d0)) :: x, y
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    REAL(kind(0.d0)), dimension(1:6) :: sigma1
    REAL(kind(0.d0)), dimension(1:3,1:3) :: stress, stress1
    
    ! Location of the point in the reference frame of the line-force
    dR(1:3) = R(1:3)-cLineCouple(1:3,n)
    x = Sum( rotLineCouple(1,1:3,n)*dR(1:3) )
    y = Sum( rotLineCouple(2,1:3,n)*dR(1:3) )

    ! Stress in the reference frame of the line-force
    vector(1:3) = 1.d0/( x + rootLineCouple(1:3,n)*y )**2
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )

    ! Stress in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotLineCouple(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotLineCouple(1:3,1:3,n) ) )
    sigma(1) = stress(1,1)
    sigma(2) = stress(2,2)
    sigma(3) = stress(3,3)
    sigma(4) = 0.5d0*( stress(2,3) + stress(3,2) )
    sigma(5) = 0.5d0*( stress(1,3) + stress(3,1) )
    sigma(6) = 0.5d0*( stress(1,2) + stress(2,1) )

  END FUNCTION LineCouple_Stress_ani

  FUNCTION LineCouple_gradStrain_ani(R, n) RESULT(gradk_Eij)
    ! Calculate strain gradient in point R(1:3) arising from the line-force couple n

    USE babel_data, ONLY : inv_CVoigt
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the line-force couple
    ! Strain gradient: gradk_Eij(i,j,k) = d[epsi(i,j)]/d[x(k)]
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3) :: gradk_Eij

    REAL(kind(0.d0)), dimension(1:3) :: dR
    REAL(kind(0.d0)) :: x, y
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    REAL(kind(0.d0)), dimension(1:6) :: sigma, sigma1, dx_epsi, dy_epsi
    REAL(kind(0.d0)), dimension(1:3,1:3) :: stress, stress1, dx_strain, dy_strain
    
    ! Location of the point in the reference frame of the line-force
    dR(1:3) = R(1:3)-cLineCouple(1:3,n)
    x = Sum( rotLineCouple(1,1:3,n)*dR(1:3) )
    y = Sum( rotLineCouple(2,1:3,n)*dR(1:3) )

    ! x derivative of the stress in the reference frame of the line-force
    vector(1:3) = -2.d0/( x + rootLineCouple(1:3,n)*y )**3
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )
    ! ... in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotLineCouple(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotLineCouple(1:3,1:3,n) ) )
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
    vector(1:3) = -2.d0*rootLineCouple(1:3,n)/( x + rootLineCouple(1:3,n)*y )**3
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )
    ! ... in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotLineCouple(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotLineCouple(1:3,1:3,n) ) )
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
    gradk_Eij(:,:,1) = inv_rotLineCouple(1,1,n)*dx_strain(:,:) &
                + inv_rotLineCouple(1,2,n)*dy_strain(:,:)
    gradk_Eij(:,:,2) = inv_rotLineCouple(2,1,n)*dx_strain(:,:) &
                + inv_rotLineCouple(2,2,n)*dy_strain(:,:)
    gradk_Eij(:,:,3) = inv_rotLineCouple(3,1,n)*dx_strain(:,:) &
                + inv_rotLineCouple(3,2,n)*dy_strain(:,:)

  END FUNCTION LineCouple_gradStrain_ani

  FUNCTION Stress_LineCouple_Interaction(sVoigt,n) RESULT(E)
    ! Interaction energy of line couple n with stress sVoigt
    USE babel_data, ONLY : inv_CVoigt
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:6), intent(in) :: sVoigt
    INTEGER, intent(in) :: n
    REAL(kind(0.d0)) :: E

    REAL(kind(0.d0)), dimension(1:6) :: eVoigt
    REAL(kind(0.d0)), dimension(1:3,1:3) :: epsi

    !   strain corresponding to this stress  + homogeneous strain
    eVoigt = MatMul( inv_CVoigt, sVoigt) 
    !   tansform to full 3x3 notation
    epsi(:,:) = Reshape( (/ eVoigt(1), 0.5d0*eVoigt(6), 0.5d0*eVoigt(5), &
                          0.5d0*eVoigt(6), eVoigt(2), 0.5d0*eVoigt(4), &
                          0.5d0*eVoigt(5), 0.5d0*eVoigt(4), eVoigt(3) /), (/3,3/) )
   !   rotate strain
   epsi(:,:) = MatMul( rotLineCouple(:,:,n), MatMul( epsi(:,:), inv_rotLineCouple(:,:,n) ) )
   !   interaction with strain of the line-force couple
   E = - epsi(1,1)*mxLineCouple(n) - epsi(2,2)*myLineCouple(n) 

  END FUNCTION Stress_LineCouple_Interaction

  SUBROUTINE RotateLineCouple(rot, out, verbose)
    ! Rotate line-force dipoles according to rotation matrix rot(1:3,1:3)
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n, i
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_rot
    
    inv_rot = Transpose( rot )

    ! Rotate definition of the line-force couples
    lLineCouple(1:3,1:nlc) = MatMul( rot(1:3,1:3), lLineCouple(1:3,1:nlc) )
    xLineCouple(1:3,1:nlc) = MatMul( rot(1:3,1:3), xLineCouple(1:3,1:nlc) )
    yLineCouple(1:3,1:nlc) = MatMul( rot(1:3,1:3), yLineCouple(1:3,1:nlc) )
    cLineCouple(1:3,1:nlc) = MatMul( rot(1:3,1:3), cLineCouple(1:3,1:nlc) )
    loop_couple1: DO n=1, nlc
            ! Beware: rotLineCouple and inv_rotLineCouple are not rank-2 tensors 
            rotLineCouple(1:3,1:3,n) = MatMul( rotLineCouple(1:3,1:3,n), inv_rot(1:3,1:3) )
            inv_rotLineCouple(1:3,1:3,n) = MatMul( rot(1:3,1:3), inv_rotLineCouple(1:3,1:3,n) )
    END DO loop_couple1

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for line-force dipoles after rotation'
            CALL PrintLineCouple(out)
    END IF

  END SUBROUTINE RotateLineCouple

  SUBROUTINE TranslateLineCouple(u, out, verbose)
    ! Add displacement u(1:3) to line-force dipole defininition
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: u
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n

    DO n=1, nlc
       cLineCouple(1:3,n) = cLineCouple(1:3,n) + u(1:3)
    END DO 

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for line-force dipoles after translation'
            CALL PrintLineCouple(out)
    END IF

  END SUBROUTINE TranslateLineCouple

END MODULE LineCouple_elasticity_ani
