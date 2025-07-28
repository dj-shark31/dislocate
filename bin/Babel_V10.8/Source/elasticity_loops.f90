MODULE Elasticity_loops

  !PRIVATE :: vecAngle

CONTAINS

  FUNCTION SolidAngle(xA, xB, xC) RESULT(Omega)
    ! Calculate solid angle (normalized by 4*pi) used to obtain the 
    ! displacement created by a dislocation triangle loop ABC at the origin
    ! (field point = origin)
    ! Ref.: Van Oosterom, A. and Strackee, J., 
    !       The Solid Angle of a Plane Triangle, 
    !       IEEE Transactions on Biomedical Engineering BME-30, 125 (1983).

    USE Math
    IMPLICIT NONE

    ! Extremities of the triangle loop
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: xA, xB, xC

    ! Solid angle (normalized by 4*pi)
    REAL(kind(0.d0)) :: omega

    REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
    REAL(kind(0.d0)), parameter :: factor=1.d0/(2.d0*pi)

    REAL(kind(0.d0)) :: rA, rB, rC, numerator, denominator

    rA = vecNorm2( xA)
    rB = vecNorm2( xB)
    rC = vecNorm2( xC)

    numerator = ScalarTripleProduct( xA, xB, xC )
    denominator = rA*rB*rC + Dot_Product( xA, xB )*rC &
        + Dot_Product( xB, xC )*rA + Dot_Product( xC, xA )*rB

    omega = factor*atan2( numerator, denominator )

  END FUNCTION SolidAngle

  FUNCTION DisloSeg_displacement_iso(xA, xB, b, nu) RESULT(u)
    ! Calculate displacement created by dislocation segment AB
    ! once the solid angle part has been removed
    ! Isotropic elastic calculation with nu Poisson coef.
    ! Ref.: Eq. (1) in Barnett, D. M. 
    !       The Displacement Field of a Triangular Dislocation Loop
    !       Philos. Mag. A, 1985, 51, 383-387 

    USE Math
    IMPLICIT NONE

    ! Extremities of the segment
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: xA, xB
    ! Burgers vector
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: b
    ! Poisson coefficient
    REAL(kind(0.d0)), intent(in) :: nu

    ! Displacement
    REAL(kind(0.d0)), dimension(1:3) :: u

    REAL(kind(0.d0)) :: rA, rB, nABnorm, tABnorm
    REAL(kind(0.d0)), dimension(1:3) ::  tAB, nAB 
    REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
    REAL(kind(0.d0)), parameter :: Zero_angle=1.d-6     ! radians

    rA = vecNorm2( xA)
    rB = vecNorm2( xB)

    ! Tangent vector
    tAB(:) = xB(:) - xA(:)
    tABnorm = vecNorm2(tAB)
    IF (tABnorm.GT.zero_angle*rA*rB) THEN
            tAB(:) = tAB(:)/tABnorm
    ELSE
            tAB(:) = 0.d0
    END IF

    ! Normal vector
    nAB(:) = CrossProduct(xA, xB)
    nABnorm = vecNorm2(nAB)
    IF (nABnorm.GT.zero_angle*rA*rB) THEN
            nAB = nAB/nABnorm
    ELSE
            nAB(:) = 0.d0
    END IF

    u(:) = ( -(1.d0-2.d0*nu)*CrossProduct(b, tAB)* &
            log( (rB + Dot_Product(xB,tAB))/(rA + Dot_Product(xA,tAB)) ) &
            + Dot_Product(b,nAB)*CrossProduct(xB/rB-xA/rA,nAB) )&
            /(8.d0*pi*(1.d0-nu))

  END FUNCTION DisloSeg_displacement_iso

  FUNCTION SolidAngle_Hirth(xA, xB, xC, r) RESULT(Omega)
    ! Calculate solid angle (normalized by 4*pi) used to obtain the 
    ! displacement created by a dislocation triangle loop ABC 
    ! Ref.: Hirth and Lothe, Eq. (5.87) - (5.90) pp. 146-147

    USE Math
    IMPLICIT NONE

    ! Extremities of the triangle loop
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: xA, xB, xC
    ! Field point
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: r

    ! Solid angle (normalized by 4*pi)
    REAL(kind(0.d0)) :: omega

    REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
    REAL(kind(0.d0)), parameter :: factor=1.d0/(4.d0*pi)

    REAL(kind(0.d0)), dimension(1:3) :: lAB, eAB, normal
    REAL(kind(0.d0)), dimension(1:3) :: rPrime, rA, rB, xP
    REAL(kind(0.d0)) :: hAB, lA2, lB2, Z, phiAB

    ! Normal of the triangle loop
    normal(:) = CrossProduct( xA(:)-xC(:), xB(:)-xC(:) )
    normal(:) = normal(:)/vecNorm2(normal(:))

    ! Line vector of the segment
    lAB(:) = xB(:) - xA(:)
    lAB(:) = lAB(:) / vecNorm2(lAB(:) )

    eAB(:) = CrossProduct(normal(:), lAB(:))

    ! Project field point on the surface
    rPrime(:) = - Dot_Product( r(:), normal(:) )*normal(:)
    xP(:) = r(:) + rPrime(:)
    rA(:) = xA(:) - r(:)
    rB(:) = xB(:) - r(:)

    hAB = Dot_Product( rPrime(:)-rA(:), eAB(:) )
    lA2 = Sum( ( rA(:) - rPrime(:) )**2 )
    lB2 = Sum( ( rB(:) - rPrime(:) )**2 )
    Z = - Dot_Product( rPrime(:), normal(:) )

    ! ==== DEBUG ===========================
    WRITE(6,'(a,3g20.12)') 'r   = ', r(1:3)
    WRITE(6,'(a,3g20.12)') 'normal   = ', normal(1:3)
    WRITE(6,'(a,3g20.12)') 'xA  = ', xA(1:3)
    WRITE(6,'(a,3g20.12)') 'xB  = ', xB(1:3)
    WRITE(6,'(a,3g20.12)') 'xP  = ', xP(1:3)
    WRITE(6,'(a,3g20.12)') 'lAB = ', lAB(1:3)
    WRITE(6,'(a,3g20.12)') 'eAB = ', eAB(1:3)
    WRITE(6,'(a,g20.12)')  'hAB = ', hAB
    WRITE(6,'(a,g20.12)')  'lA2 = ', lA2
    WRITE(6,'(a,g20.12)')  'lB2 = ', lB2
    WRITE(6,'(a,g20.12)')  'Z = ', Z
    WRITE(6,'(a,g20.12)')  'lA2-hAB**2 = ', lA2-hAB**2
    WRITE(6,'(a,g20.12)')  'lB2-hAB**2 = ', lB2-hAB**2
    WRITE(6,*)
    ! ==== DEBUG ===========================

    phiAB = vecAngle(xA-xP, xB-xP, normal)
    IF (phiAB.LT.0.d0) phiAB=phiAB+2.d0*pi

    Omega = factor*( ( phiAB-pi )*sign(1.d0,Z) & 
        + atan( hAB/Z*Sqrt( (lA2+Z**2)/(lA2-hAB**2) ) ) &
        + atan( hAB/Z*Sqrt( (lB2+Z**2)/(lB2-hAB**2) ) ) )


  END FUNCTION SolidAngle_Hirth

  FUNCTION vecAngle(u1, u2, n) RESULT(phi)
    ! Angles between vectors u1 and u2, oriented according to normal n
    
    USE Math
    IMPLICIT NONE

    REAL(kind(0.d0)), dimension(1:3), intent(in) :: u1, u2, n
    REAL(kind(0.d0)) :: phi

    !! ==== DEBUG ===========================
    !WRITE(6,'(a,3g20.12)') 'u1 = ', u1(1:3)
    !WRITE(6,'(a,3g20.12)') 'u2 = ', u2(1:3)
    !WRITE(6,'(a,g20.12)') 'vecNorm2(u1) = ', vecNorm2(u1)
    !WRITE(6,'(a,g20.12)') 'vecNorm2(u2) = ', vecNorm2(u2)
    !WRITE(6,*)
    !! ==== DEBUG ===========================

    phi = sign( acos( dot_product(u1, u2)/ ( vecnorm2(u1)*vecnorm2(u2) ) ), ScalarTripleProduct(u1,u2,n) )

  END FUNCTION vecAngle

END MODULE Elasticity_loops

