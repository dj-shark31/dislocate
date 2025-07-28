MODULE StrainFromDisplacementModule

  ! Euler Lagrange definition used for the strain when EulerLagrange=.true.
  ! Green Lagrange definition used otherwise (default)
  LOGICAL, save :: EulerLagrange=.FALSE.

  PRIVATE :: Mat3Inv

CONTAINS

  SUBROUTINE StrainFromDisplacement(im, gradDisplacement, eAtomicStrain, out)
     ! Calculate strain from displacement gradient
     ! Green-Lagrange or Euler-Lagrange definition depending on logical EulerLagrange

     USE Babel_data, ONLY : verbosity
     IMPLICIT NONE

     ! Number of atoms
     INTEGER, intent(in) :: im
     ! Displacement gradient of each atom
     REAL(kind(0.d0)), dimension(:,:,:), intent(in) :: gradDisplacement
     ! Strain for each atom (Voigt notation)
     REAL(kind(0.d0)), dimension(:,:), intent(out) :: eAtomicStrain
     ! Output unit 
     INTEGER, intent(in), optional :: out

     INTEGER :: i
     REAL(kind(0.d0)), dimension(1:3,1:3) :: epsi, dudx, Tdudx, Amat, Bmat

     DO i=1, im
        dudx(:,:) = gradDisplacement(1:3,1:3,i)
        Tdudx(:,:) = Transpose(dudx(:,:))
        epsi(:,:) = 0.5d0*( dudx(:,:) + Tdudx(:,:) + MatMul( Tdudx(:,:), dudx(:,:) ) )
        IF (EulerLagrange) THEN
                ! Euler-Lagrange definition is used for the strain
                Amat(:,:) = dudx(:,:) 
                Amat(1,1) = Amat(1,1) + 1.d0
                Amat(2,2) = Amat(2,2) + 1.d0
                Amat(3,3) = Amat(3,3) + 1.d0
                CALL Mat3Inv(Amat, Bmat)
                epsi = MatMul( Transpose(Bmat), MatMul( epsi, Bmat ) )
        END IF
        eAtomicStrain(1,i) = epsi(1,1)
        eAtomicStrain(2,i) = epsi(2,2)
        eAtomicStrain(3,i) = epsi(3,3)
        eAtomicStrain(4,i) = epsi(3,2) + epsi(2,3)
        eAtomicStrain(5,i) = epsi(3,1) + epsi(1,3)
        eAtomicStrain(6,i) = epsi(2,1) + epsi(1,2)
     END DO

  END SUBROUTINE StrainFromDisplacement

  SUBROUTINE Mat3Inv(A, B)
    ! return matrix B which is the inverse of matrix A
    ! A and B are 3x3 matrices

    implicit none

    REAL(kind(0.d0)), dimension(3,3), intent(in) :: A
    REAL(kind(0.d0)), dimension(3,3), intent(out) :: B

    REAL(kind(0.d0)) :: invdet

       
    b(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
    b(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
    b(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
       
    b(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
    b(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
    b(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)
       
    b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
    b(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    invdet = 1.d0/( a(1,1)*b(1,1) + a(1,2)*b(2,1) + a(1,3)*b(3,1) )
    b(:,:)=b(:,:)*invdet

  END SUBROUTINE Mat3Inv

END MODULE StrainFromDisplacementModule
