MODULE dDipoleModule

  ! Definition of the dislocation dipole
  LOGICAL, save :: l_DDipole      ! True if dislo dipoles defined
  INTEGER, save :: nDDipole       ! Number of dislocation dipoles in the simulation box
  INTEGER, parameter :: max_nDDipole=200    ! Maximal value for nDDipole
  REAL(kind(0.d0)), dimension(1:3,1:max_nDDipole), save :: lDDipole    ! Line vector
  REAL(kind(0.d0)), dimension(1:3,1:max_nDDipole), save :: bDDipole    ! Burgers vector (in A)
  REAL(kind(0.d0)), dimension(1:3,1:max_nDDipole), save :: c1DDipole, c2DDipole      ! Point belonging to the dislo +b and -b (in A)
  REAL(kind(0.d0)), dimension(1:3,1:3,1:max_nDDipole), save :: rotDDipole, inv_rotDDipole ! Matrix for changing ref. frame

  ! Variable set to true if only one dislocation system (one line direction) is
  ! present => allow energy calculation
  LOGICAL, save :: one_DDipole_system

  !---------------------------------------------------------
  ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
  !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
  !   and the equation det[a(i;j)]=0 has 6 solutions for the variable p: root(1:6) 
  COMPLEX(kind(0.d0)), dimension(1:6,1:max_nDDipole), save, private :: root
  ! Matrices used to calculate displacement and stress field created by a dislo
  COMPLEX(kind(0.d0)), dimension(1:3,1:6,1:max_nDDipole), save, private :: displacement_matrix
  COMPLEX(kind(0.d0)), dimension(1:6,1:6,1:max_nDDipole), save, private :: stress_matrix

  ! Angular dependence of the energy
  REAL(kind(0.d0)), dimension(1:max_nDDipole), save, public :: DDipole_angular_energy
  ! Prelogarithmic energy factor used to calculate dislocation elastic energy
  REAL(kind(0.d0)), dimension(1:max_nDDipole), save, public :: DDipole_prelog_energy_factor
  ! Matrix used to calculate interaction energy between dislocations
  COMPLEX(kind(0.d0)), dimension(1:3,1:6,1:max_nDDipole), save, private :: interaction

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE InitDDipole()

    IMPLICIT NONE

    l_DDipole = .FALSE.
    nDDipole = 0
    lDDipole(:,:) = 0.d0
    bDDipole(:,:) = 0.d0
    c1DDipole(:,:) = 0.d0 
    c2DDipole(:,:) = 0.d0 
    rotDDipole(:,:,:) = 0.d0
    inv_rotDDipole(:,:,:) = 0.d0

  END SUBROUTINE InitDDipole

  SUBROUTINE ReadDDipole(inp)

    IMPLICIT NONE
    INTEGER, intent(in) :: inp

    NAMELIST /DDipole/ nDDipole, lDDipole, bDDipole, c1DDipole, c2DDipole

    CALL InitDDipole()

    READ(inp,nml=ddipole)
    l_ddipole=.true.

    ! Check number of dislo dipole
    IF (nDDipole.LE.0) THEN
            CALL InitDDipole()
            RETURN
    ELSEIF (nDDipole.GT.max_nDDipole) THEN
            WRITE(0,'(a,i0)') 'Maximal number of dislocation dipoles (max_nDDipole): ', max_nDDipole
            WRITE(0,'(a,i0)') 'Current number of dislocation dipoles (nDDipole): ', nDDipole
            WRITE(0,'(a)') 'you need to change value in disloDipole.f90 and to re-compile'
            STOP '< ReadDDipole >'
    END IF

  END SUBROUTINE ReadDDipole
  
  SUBROUTINE PrintDDipole(out)
    ! Print dislocation dipole definition on output unit out

    IMPLICIT NONE
    INTEGER, intent(in) :: out

    INTEGER :: n, i

    IF (.NOT.l_DDipole) THEN
            WRITE(out,'(a)') ' No dislocation dipole defined'
            WRITE(out,*)
            RETURN
    END IF

    DO n=1, nDDipole
            WRITE(out,'(a,i0)') '  Dislocation dipole ', n
            WRITE(out,'(a,3g14.6)') '    line vector:         ', lDDipole(1:3,n)
            WRITE(out,'(a,3g14.6)') '    Burgers vector (A):  ', bDDipole(1:3,n)
            WRITE(out,'(a,3g14.6)') '    center 1 (A) for +b: ', c1DDipole(1:3,n)
            WRITE(out,'(a,3g14.6)') '    center 2 (A) for -b: ', c2DDipole(1:3,n)
            WRITE(out,'(a)') '    corresponding rotation matrix:'
            DO i=1, 3
               WRITE(out,'(6x,3f18.4)') rotDDipole(i,1:3,n)
            END DO
    END DO
    WRITE(out,*)

  END SUBROUTINE PrintDDipole

  SUBROUTINE WriteDDipole(out)
    ! Print dislocation dipole definition on output unit out

    IMPLICIT NONE
    INTEGER, intent(in) :: out

    INTEGER :: n

    IF (.NOT.l_DDipole) RETURN
    WRITE(out,'(a)') "#  1: dislocation dipole number n"
    WRITE(out,'(a)') "#  2: dipole Burgers vector, bx"
    WRITE(out,'(a)') "#  3:                        by"
    WRITE(out,'(a)') "#  4:                        bz"
    WRITE(out,'(a)') "#  5: dipole center 1, x"
    WRITE(out,'(a)') "#  6:                  y"
    WRITE(out,'(a)') "#  7:                  z"
    WRITE(out,'(a)') "#  8: dipole center 2, x"
    WRITE(out,'(a)') "#  9:                  y"
    WRITE(out,'(a)') "# 10:                  z"
    DO n=1, nDDipole
            WRITE(out,'(i0,1x,9g14.6)') n, bDDipole(1:3,n), c1DDipole(1:3,n), c2DDipole(1:3,n)
    END DO
    WRITE(out,*)

  END SUBROUTINE WriteDDipole

  SUBROUTINE ScaleDDipole(a)
    ! Multiply distances by a for dislocation dipole definition

    IMPLICIT NONE
    REAL(kind(0.d0)), intent(in) :: a

    ! Multiply dislocation dipole Burgers vector and centers
    bDDipole(:,1:nDDipole)  = a*bDDipole(:,1:nDDipole)
    c1DDipole(:,1:nDDipole) = a*c1DDipole(:,1:nDDipole)
    c2DDipole(:,1:nDDipole) = a*c2DDipole(:,1:nDDipole)

  END SUBROUTINE ScaleDDipole

  SUBROUTINE RearrangeDDipole()

    USE babel_data, ONLY : xImages, yImages, zImages, at, &
        distance_zero, matid
    USE math
    IMPLICIT NONE

    INTEGER :: n
    REAL(kind(0.d0)) :: inv_norm, scalar
    REAL(kind(0.d0)), dimension(1:3) :: cut, ex, ey, ez

    IF (nDDipole.LE.0) RETURN

    DO n=1, nDDipole
 
       ! Normalize dislocation line vectors
       inv_norm = 1.d0/sqrt( Sum( lDDipole(1:3,n)**2 ) )
       lDDipole(1:3,n) = inv_norm*lDDipole(1:3,n)

       ! Check if all dipole have the same direction and change the sign if necessary
       one_DDipole_system=.true.
       IF (n.NE.1) THEN 
               scalar = Sum( lDDipole(1:3,1)*lDDipole(1:3,n) )
               IF (Abs(scalar+1.d0).LE.distance_zero) THEN 
                       ! Dipole n belong to the same system as dipole 1 
                       ! but we need to invert its orientation
                       lDDipole(1:3,n) = - lDDipole(1:3,n)
                       bDDipole(1:3,n) = - bDDipole(1:3,n)
               ELSE IF (Abs(scalar-1.d0).GT.distance_zero) THEN
                       one_DDipole_system=.FALSE.
               END IF
       END IF
       
       ! Define cut direction
       cut(1:3) = c2DDipole(1:3,n) - c1DDipole(1:3,n)

       ! Make it orthogonal to line direction
       scalar = Sum( lDDipole(:,n)*cut(:) )
       cut(:) = cut(:) - scalar*lDDipole(:,n)

       ! Normalize it and check that it is not parallel to line direction
       scalar = Sqrt( Sum( cut(:)**2 ) )
       IF (scalar.LE.distance_zero) THEN
               WRITE(0,'(a,i0)') "Dislocation dipole ", n
               WRITE(0,'(a)') "The vector joining dislocation centers is parallel to line direction"
               STOP '< RearrangeDDipole >'
       END IF
       cut(:) = cut(:)/scalar

       ! Orientate the crystal
       !  z: dislo vector line
       !  y: glide plane normal
       !  x: cutting direction corresponding to the discontinuity
       ez(1:3) = lDDipole(1:3,n)
       ex(1:3) = cut(:)
       ey(1:3) = CrossProduct( ez(1:3), ex(1:3) ) 

       ! Matrix to change the orientation of the crystal from cartesian coordinates
       ! (cubic axes) to the ref. frame of the dislocation
       rotDDipole(1,1:3,n) = ex(1:3)
       rotDDipole(2,1:3,n) = ey(1:3)
       rotDDipole(3,1:3,n) = ez(1:3)

       ! Inverse matrix
       inv_rotDDipole(1:3,1:3,n) = Transpose( rotDDipole(1:3,1:3,n) )
       
       ! Check that rotDislo is a rotation matrix
       scalar = MatNorm2( MatMul( inv_rotDDipole(:,:,n), rotDDipole(:,:,n) ) - matId )
       IF (scalar.GT.1.d-6) THEN
               WRITE(0,'(a,i0,a)') 'Rotation matrix for dislo dipole ', n, ' is not a unitary matrix'
               STOP '< RearrangeDDipole >'
       END IF
    END DO

    !------------------------------------
    ! Check that periodic boundary conditions are compatible with 
    ! line directions and Burgers circuits for all systems
    IF (xImages.OR.yImages.OR.zImages) THEN
            DO n=1, nDDipole ! Loop on all systems
               IF ( ( xImages .AND. ( Colinear(at(:,1),lDDipole(:,n)) ) ).OR. &
                    ( yImages .AND. ( Colinear(at(:,2),lDDipole(:,n)) ) ).OR. &
                    ( zImages .AND. ( Colinear(at(:,3),lDDipole(:,n)) ) ) ) THEN
                       WRITE(0,'(a,i0,a)') 'Line direction of dislocations dipole ', &
                                n, ' is colinear with a periodicity vector'
                       STOP '< RearrangeDDipole >'
               END IF
            END DO
    END IF
  END SUBROUTINE RearrangeDDipole

  SUBROUTINE RotateDDipole(rot, out, verbose)
    ! Rotate dislocation dipoles according to rotation matrix rot(1:3,1:3)
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_rot
    
    inv_rot = Transpose( rot )

     lDDipole(:,1:nDDipole) = MatMul( rot(:,:),  lDDipole(:,1:nDDipole) )
     bDDipole(:,1:nDDipole) = MatMul( rot(:,:),  bDDipole(:,1:nDDipole) )
    c1DDipole(:,1:nDDipole) = MatMul( rot(:,:), c1DDipole(:,1:nDDipole) )
    c2DDipole(:,1:nDDipole) = MatMul( rot(:,:), c2DDipole(:,1:nDDipole) )
    DO n=1, nDDipole
            ! Beware: rotDDipole and inv_rotDDipole are not rank-2 tensors 
            rotDDipole(:,:,n) = MatMul( rot(:,:), MatMul( rotDDipole(:,:,n), inv_rot(:,:) ) )
            inv_rotDDipole(:,:,n) = MatMul( rot(:,:), inv_rotDDipole(:,:,n) )
    END DO

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for dislocation dipoles after rotation'
            CALL PrintDDipole(out)
    END IF

  END SUBROUTINE RotateDDipole

  SUBROUTINE TranslateDDipole(u, out, verbose)
    ! Add displacement u(1:3) to dislocation dipole defininition
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: u
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n

    DO n=1, nDDipole
       c1DDipole(1:3,n) = c1DDipole(1:3,n) + u(1:3)
       c2DDipole(1:3,n) = c2DDipole(1:3,n) + u(1:3)
    END DO 

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for dislocation dipoles after translation'
            CALL PrintDDipole(out)
    END IF

  END SUBROUTINE TranslateDDipole

  SUBROUTINE InitDDipole_ani(out)
    ! Initialization of all quantities used for anisotropic elastic calculation     
    !   for a dislocation dipole

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
    REAL(kind(0.d0)) :: norm, d
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt_noisy
    REAL(kind(0.d0)), dimension(1:3) :: dR

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
            STOP '< InitDDipole_ani >'
    END IF

    ! Initialization
    elastic_C = Unpack_CVoigt(CVoigt_noisy)

    DO n=1, nDDipole

            ! Rotate elastic constants in the dislocation axes
            elastic_Cp = Rotate_Elastic_Constants(elastic_C, rotDDipole(1:3,1:3,n) )
            ! Rotate Burgers vector in dislo axes: burgers(1) = edge component (if non 0)
            !                                      burgers(2) = 0
            !                                      burgers(3) = screw component
            burgers(1:3) = MatMul( rotDDipole(1:3,1:3,n),bDDipole(1:3,n) )
            
            ! No line-force for a dislocation
            force(1:3) = 0.d0

            IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
                    WRITE(out,'(a,i0)') ' Initialization for dislocation dipole ', n
                    WRITE(out,'(a,3g14.6)') '    Burgers vector in rotated axes (A): ', &
                        burgers(1:3)
                    WRITE(out,'(a,3g14.6)') '    line direction in rotated axes : ', &
                        MatMul( rotDDipole(1:3,1:3,n), lDDipole(1:3,n) )
                    WRITE(out,'(a,3g14.6)') '    cut direction in rotated axes : ', &
                        MatMul( rotDDipole(1:3,1:3,n), c2DDipole(1:3,n)-c1DDipole(1:3,n) )
                    WRITE(out,'(a)') '  elastic constants in rotated axes'
                    CALL Print_Elastic_Constants(elastic_Cp,out)
            END IF


            IF (Present(out)) THEN
                    CALL InitStroh(elastic_Cp, root(:,n), elastic_A, elastic_B, elastic_L, out)
                    CALL Build_DStroh(burgers, force, root(:,n), elastic_A, elastic_B, elastic_L, &
                        displacement_matrix(:,:,n), stress_matrix(:,:,n), &
                        DDipole_angular_energy(n), DDipole_prelog_energy_factor(n), interaction(:,:,n), out)
            ELSE
                    CALL InitStroh(elastic_Cp, root(:,n), elastic_A, elastic_B, elastic_L)
                    CALL Build_DStroh(burgers, force, root(:,n), elastic_A, elastic_B, elastic_L, &
                        displacement_matrix(:,:,n), stress_matrix(:,:,n), &
                        DDipole_angular_energy(n), DDipole_prelog_energy_factor(n), interaction(:,:,n))
            END IF

            ! Displacement in cartesian coordinate
            displacement_matrix(:,:,n) = MatMul( inv_rotDDipole(:,:,n), displacement_matrix(:,:,n) )

            IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
                    b2 = Sum(bDDipole(1:3,n)**2)

                    ! prelogarithmic energy factor: E = K*ln(R/r)
                    WRITE(out,*)
                    WRITE(out,'(a)') '  Elastic energy of the isolated dipole: Edipole = 2*E0 + K0 ln(|d|/rc)'
                    WRITE(out,'(a)') '  prelogarithmic energy factor:'
                    WRITE(out,'(a,g20.12,a)') '     K0 = ', 2.d0*DDipole_prelog_energy_factor(n), &
                        '  -  units = elastic constant * distance^2 (GPa.A^2)'
                    WRITE(out,'(a,g20.12,a)') '        = ', 2.d0*factorE*DDipole_prelog_energy_factor(n), &
                        '  -  units = energy / distance (eV.A^-1)'
                    WRITE(out,'(a,g20.12,a)') '        = ', 2.d0*DDipole_prelog_energy_factor(n)/b2*2*pi, &
                        '*b^2/2pi  -  units = elastic constant (GPa)'

                    ! angular dependence of the energy
                    WRITE(out,'(a)') '  angular dependence of the energy'
                    WRITE(out,'(a,g20.12,a)') '    E0 = ', DDipole_angular_energy(n), &
                        '  -  units = elastic constant * distance^2 (GPa.A^2)'
                    WRITE(out,'(a,g20.12,a)') '       = ', factorE*DDipole_angular_energy(n), &
                        '  -  units = energy / distance (eV.A^-1)'
                    WRITE(out,'(a,g20.12,a)') '       = ', DDipole_angular_energy(n)/b2*2*pi, &
                        '*b^2/2pi  -  units = elastic constant (GPa)'

                    ! separation distance of the dipole
                    !d = sqrt( Sum( ( c1DDipole(:,n)-c2DDipole(:,n) )**2 ) )
                    dR(:) = c1DDipole(:,n) - c2DDipole(:,n) &
                            - Sum( ( c1DDipole(:,n) - c2DDipole(:,n) )*lDDipole(:,n) )*lDDipole(:,n)
                    d = Sqrt( Sum( dR(:)**2 ) )
                    WRITE(out,'(a)') '  separation distance of the dipole'
                    WRITE(out,'(a,g20.12,a)') '    d = ', d,  '  -  units = distance (A)'

                    ! Prelogarithmic tensor appearing in elastic energy: E = 1/2 bi Kij bj ln(R/rc)
                    DO i=1, 3 ; DO j=1, 3
                        Kij(i,j) = aImag( Sum( elastic_L(i,1:3)*elastic_L(j,1:3) ) &
                                - Sum( elastic_L(i,4:6)*elastic_L(j,4:6) ) )
                    END DO ; END DO
                    Kij(:,:) = Kij(:,:)/(2.d0*pi)
                    Kij(:,:) = MatMul( Transpose( rotDDipole(:,:,n) ), MatMul( &
                        Kij(:,:), rotDDipole(:,:,n) ) )
                    WRITE(out,'(a)') '  Stroh prelogarithmic tensor appearing in interaction energy: E = bi^1 Kij bj^2 ln(R/r)  (original axes)'
                    WRITE(out,'(a,3g20.12,a)') '           | ', Kij(1,1:3), ' | '
                    WRITE(out,'(a,3g20.12,a)') '     Kij = | ', Kij(2,1:3), ' | &
                        & units = elastic constant (GPa)'
                    WRITE(out,'(a,3g20.12,a)') '           | ', Kij(3,1:3), ' | '
                    K0 = factorE*Sum( bDDipole(:,n) * MatMul( Kij(:,:), bDDipole(:,n) ) )
                    WRITE(out,'(a,g20.12,a)') '    check for dipole self-energy: K0 ?= bi Kij bj = ', K0, &
                        '  -  units = energy / distance (eV.A^-1)'
                    WRITE(out,*)

            END IF

    END DO
            
  END SUBROUTINE InitDDipole_ani

  FUNCTION DDipole_Displacement_ani(R, n) RESULT(u)
    ! Calculate total displacement in point R(1:3) arising from the dislocation dipole n
    !   ref.: Hirth and Lothe, p.445, Eq.13.91

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the dislocation dipole
    REAL(kind(0.d0)), dimension(1:3) :: u

    REAL(kind(0.d0)), dimension(2) :: ds1, ds2
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    
    ! Location of the point in the reference frame of the dislo
    ds1(:) = MatMul( rotDDipole(1:2,:,n), R(:)-c1DDipole(:,n))
    ds2(:) = MatMul( rotDDipole(1:2,:,n), R(:)-c2DDipole(:,n))

    !WRITE(6,*) 'ds1: ', ds1(:)
    !WRITE(6,*) 'ds2: ', ds2(:)

    IF (Sum(ds1(:)**2).LT.1.d-8) THEN
            WRITE(0,'(a,3g14.6)') 'Point R(1:3) = ', R(1:3)
            WRITE(0,'(a,3g14.6)') 'Dislocation center c1DDipole(1:3) = ', c1DDipole(1:3,n)
            WRITE(0,'(2(a,g14.6))') ' Distance from dislocation line: x = ', ds1(1), '  - y = ', ds1(2)
            WRITE(0,'(a)') '  => cannot calculate displacement'
            !STOP "< DDipole_Displacement_ani >"
            WRITE(0,*)
            u(:)=0.d0
            RETURN
    ELSE IF (Sum(ds2(:)**2).LT.1.d-8) THEN
            WRITE(0,'(a,3g14.6)') 'Point R(1:3) = ', R(1:3)
            WRITE(0,'(a,3g14.6)') 'Dislocation center c2DDipole(1:3) = ', c2DDipole(1:3,n)
            WRITE(0,'(2(a,g14.6))') ' Distance from dislocation line: x = ', ds2(1), '  - y = ', ds2(2)
            WRITE(0,'(a)') '  => cannot calculate displacement'
            !STOP "< DDipole_Displacement_ani >"
            WRITE(0,*)
            u(:)=0.d0
    END IF

    ! Displacement in the reference frame of the dislo
    vector(1:3) = log( ds1(1) + root(1:3,n)*ds1(2) ) - log( ds2(1) + root(1:3,n)*ds2(2) )
    u(:) = Dble( MatMul(displacement_matrix(:,1:3,n), vector(:)) )

  END FUNCTION DDipole_Displacement_ani

  FUNCTION DDipole_Stress_ani(R, n, onlyDislo) RESULT(sigma)
    ! Calculate stress in point R(1:3) arising from the dislocation dipole n
    !   if onlyDislo is present and equal to 1 or 2, then only dislocation 1 or 2 is considered 
    !   ref.: Hirth and Lothe, p.445, Eq.13.92

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the dislocation
    INTEGER, intent(in), optional :: onlyDislo  ! only one dislo (1 or 2) is considered
    ! Stress in Voigt notation
    REAL(kind(0.d0)), dimension(1:6) :: sigma

    REAL(kind(0.d0)), dimension(1:3) :: dR1, dR2
    REAL(kind(0.d0)) :: x1, y1, x2, y2
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    REAL(kind(0.d0)), dimension(1:6) :: sigma1
    REAL(kind(0.d0)), dimension(1:3,1:3) :: stress, stress1
    
    ! Location of the point in the reference frame of the dislo
    dR1(1:3) = R(1:3)-c1DDipole(1:3,n)
    x1 = Sum( rotDDipole(1,1:3,n)*dR1(1:3) )
    y1 = Sum( rotDDipole(2,1:3,n)*dR1(1:3) )
    dR2(1:3) = R(1:3)-c2DDipole(1:3,n)
    x2 = Sum( rotDDipole(1,1:3,n)*dR2(1:3) )
    y2 = Sum( rotDDipole(2,1:3,n)*dR2(1:3) )

    ! Stress in the reference frame of the dislo
    IF (Present(onlyDislo)) THEN
            IF (onlyDislo.EQ.1) THEN
                    vector(1:3) =  1.d0/(x1+root(1:3,n)*y1)
            ELSEIF (onlyDislo.EQ.2) THEN
                    vector(1:3) = -1.d0/(x2+root(1:3,n)*y2)
            ELSE
                    vector(1:3) =  1.d0/(x1+root(1:3,n)*y1) - 1.d0/(x2+root(1:3,n)*y2)
            END IF
    ELSE
            vector(1:3) =  1.d0/(x1+root(1:3,n)*y1) - 1.d0/(x2+root(1:3,n)*y2)
    END IF
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )

    ! Stress in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotDDipole(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotDDipole(1:3,1:3,n) ) )
    sigma(1) = stress(1,1)
    sigma(2) = stress(2,2)
    sigma(3) = stress(3,3)
    sigma(4) = 0.5d0*( stress(2,3) + stress(3,2) )
    sigma(5) = 0.5d0*( stress(1,3) + stress(3,1) )
    sigma(6) = 0.5d0*( stress(1,2) + stress(2,1) )

  END FUNCTION DDipole_Stress_ani

  FUNCTION DDipole_gradStrain_ani(R, n) RESULT(gradk_Eij)
    ! Calculate stain gradient in point R(1:3) arising from the dislocation dipole n

    USE babel_data, ONLY : inv_CVoigt    
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: R        ! (A)
    INTEGER, intent(in) :: n    ! Index of the dislocation
    ! Strain gradient: gradk_Eij(i,j,k) = d[epsi(i,j)]/d[x(k)]    
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3) :: gradk_Eij

    REAL(kind(0.d0)), dimension(1:3) :: dR1, dR2
    REAL(kind(0.d0)) :: x1, y1, x2, y2
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector
    REAL(kind(0.d0)), dimension(1:6) :: sigma, sigma1, dx_epsi, dy_epsi
    REAL(kind(0.d0)), dimension(1:3,1:3) :: stress, stress1, dx_strain, dy_strain
    
    ! Location of the point in the reference frame of the dislo
    dR1(1:3) = R(1:3)-c1DDipole(1:3,n)
    x1 = Sum( rotDDipole(1,1:3,n)*dR1(1:3) )
    y1 = Sum( rotDDipole(2,1:3,n)*dR1(1:3) )
    dR2(1:3) = R(1:3)-c2DDipole(1:3,n)
    x2 = Sum( rotDDipole(1,1:3,n)*dR2(1:3) )
    y2 = Sum( rotDDipole(2,1:3,n)*dR2(1:3) )

    ! x derivative of the stress in the reference frame of the dislo
    vector(1:3) =  -1.d0/(x1+root(1:3,n)*y1)**2 + 1.d0/(x2+root(1:3,n)*y2)**2
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )
    ! ... in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotDDipole(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotDDipole(1:3,1:3,n) ) )
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
    vector(1:3) =  -root(1:3,n)/(x1+root(1:3,n)*y1)**2 + root(1:3,n)/(x2+root(1:3,n)*y2)**2
    sigma1(1:6) = Dble( MatMul(stress_matrix(1:6,1:3,n),vector(1:3)) )
    stress1 = Reshape( (/ sigma1(1), sigma1(6), sigma1(5), &
                          sigma1(6), sigma1(2), sigma1(4), &
                          sigma1(5), sigma1(4), sigma1(3) /), (/3, 3/) )
    ! ... in cartesian coordinates
    stress(1:3,1:3) = MatMul( inv_rotDDipole(1:3,1:3,n), &
        MatMul( stress1(1:3,1:3), rotDDipole(1:3,1:3,n) ) )
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
    gradk_Eij(:,:,1) = inv_rotDDipole(1,1,n)*dx_strain(:,:) &
                + inv_rotDDipole(1,2,n)*dy_strain(:,:)
    gradk_Eij(:,:,2) = inv_rotDDipole(2,1,n)*dx_strain(:,:) &
                + inv_rotDDipole(2,2,n)*dy_strain(:,:)
    gradk_Eij(:,:,3) = inv_rotDDipole(3,1,n)*dx_strain(:,:) &
                + inv_rotDDipole(3,2,n)*dy_strain(:,:)

  END FUNCTION DDipole_gradStrain_ani

  FUNCTION DDipole_Self_Energy(n, rc) RESULT(E)
    ! Calculate self elastic energy of the dislocation dipole n
    ! with the core cutoff rc

    IMPLICIT NONE
    INTEGER, intent(in) :: n 
    REAL(kind(0.d0)), intent(in) :: rc
    REAL(kind(0.d0)) :: E
    
    REAL(kind(0.d0)) :: d
    REAL(kind(0.d0)), dimension(1:3) :: dR

    ! separation distance of the dipole
    !d = sqrt( Sum( ( c1DDipole(:,n)-c2DDipole(:,n) )**2 ) )
    dR(:) = c1DDipole(:,n) - c2DDipole(:,n) &
            - Sum( ( c1DDipole(:,n) - c2DDipole(:,n) )*lDDipole(:,n) )*lDDipole(:,n)
    d = Sqrt( Sum( dR(:)**2 ) )

    ! Interaction energy
    E = 2.d0*DDipole_prelog_energy_factor(n)*log(d/rc) &
        + 2.d0* DDipole_angular_energy(n)

  END FUNCTION DDipole_Self_Energy

  FUNCTION DDipole_Interaction_Energy(n, m, shift) RESULT(E)
    ! Calculate interaction elastic energy between dislocation dipoles n and m
    !   dislocations have to be parallel

    IMPLICIT NONE
    INTEGER, intent(in) :: n, m
    REAL(kind(0.d0)), dimension(1:3), intent(in), optional :: shift
    REAL(kind(0.d0)) :: E
    
    REAL(kind(0.d0)), dimension(1:3) :: dR, Kn, bm
    REAL(kind(0.d0)) :: xp, yp, xq, yq, Ax, Ay
    COMPLEX(kind(0.d0)), dimension(1:3) :: vector

    ! Position of dislo p=(m,1) in the reference frame of dislo (n,1)
    dR(:) = c1DDipole(:,m) - c1DDipole(:,n)
    IF (Present(shift)) dR(:) = dR(:) + shift(:)
    xp = Sum( rotDDipole(1,:,n)*dR(:) )
    yp = Sum( rotDDipole(2,:,n)*dR(:) )

    ! Position of dislo q=(m,2) in the reference frame of dislo (n,1)
    dR(:) = c2DDipole(:,m) - c1DDipole(:,n)
    IF (Present(shift)) dR(:) = dR(:) + shift(:)
    xq = Sum( rotDDipole(1,:,n)*dR(:) )
    yq = Sum( rotDDipole(2,:,n)*dR(:) )

    ! Position of dislo m=(n,2) in the reference frame of dislo (n,1)
    dR(:) = c2DDipole(:,n) - c1DDipole(:,n)
    Ax = Sum( rotDDipole(1,:,n)*dR(:) )
    Ay = Sum( rotDDipole(2,:,n)*dR(:) )

    ! Rotate Burgers vector of dipole m in the reference frame of dipole n
    bm(:) = MatMul( rotDDipole(:,:,n),bDDipole(:,m) )

    ! Interaction energy
    vector(1:3) = log( ( 1.d0 - ( Ax + root(1:3,n)*Ay )/( xp + root(1:3,n)*yp ) ) &
        / ( 1.d0 - ( Ax + root(1:3,n)*Ay )/( xq + root(1:3,n)*yq ) ) )
    Kn(:) = Dble(MatMul(interaction(:,1:3,n),vector(1:3)))
    E = Sum( Kn(:)*bm(:) )

  END FUNCTION DDipole_Interaction_Energy

END MODULE dDipoleModule
