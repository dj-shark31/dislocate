MODULE NyeTensorModule

   IMPLICIT NONE

   ! Use periodic boundary conditions to calculate Nye tensor
   LOGICAL :: perioNyeTensor

   ! Nye tensor for each atom
   REAL(kind(0.d0)), dimension(:,:,:), allocatable :: nyeAtomic

  INTEGER, parameter, private :: verbosity_max=3

  PRIVATE :: mat3inv

CONTAINS

  SUBROUTINE InitNyeTensor(imm)

     IMPLICIT NONE

     ! Maximal number of atoms
     INTEGER, intent(in) :: imm

     perioNyeTensor=.true.

     IF (Allocated(nyeAtomic)) Deallocate(nyeAtomic)
     ALLOCATE(nyeAtomic(1:3,1:3,1:imm))
     nyeAtomic(:,:,:) = 0.d0

  END SUBROUTINE InitNyeTensor

  SUBROUTINE BuildNyeTensor(im, xp0, du, at0,  nNeigh, iNeigh, out)
     ! Calculate Nye tensor from finite differences using
     ! displacement gradients of neighbours for each atom

     USE Babel_data, ONLY : verbosity
     IMPLICIT NONE

     ! Number of atoms
     INTEGER, intent(in) :: im
     ! Atoms coordinates in reference configuration
     REAL(kind(0.d0)), dimension(:,:), intent(in) :: xp0
     ! Displacement gradient of each atom
     !   du(i,j, n) = d U_i / d X_j for atom n
     REAL(kind(0.d0)), dimension(:,:,:), intent(in) :: du
     ! Periodicity vectors of the reference structure: at0(1:3,1), at0(1:3,2), ...
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at0
     ! Number of neighbours for each atom 
     INTEGER, dimension(:), intent(in) :: nNeigh
     ! Neighbour indexes for each atom i: iNeigh(1:nNeigh(i), i)
     INTEGER, dimension(:,:), intent(in) :: iNeigh
     ! Output unit 
     INTEGER, intent(in), optional :: out

     REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at0
     REAL(kind(0.d0)), dimension(:,:), allocatable :: R0
     REAL(kind(0.d0)), dimension(:,:,:), allocatable :: dG
     REAL(kind(0.d0)), dimension(1:3) :: ds0
     REAL(kind(0.d0)), dimension(1:3,1:3) :: Amat, Bmat, invBmat
     REAL(kind(0.d0)), dimension(1:3,1:3,1:3) :: Vmat, dGdx
     REAL(kind(0.d0)), dimension(1:3,1:3,1:im) :: g
     INTEGER :: i, j, l, n, nMax

     ! Inverse periodicity vectors
     IF (perioNyeTensor) THEN
             CALL Mat3Inv(at0, inv_at0)
     ELSE
             inv_at0(:,:)=0.d0
     END IF

     ! Calculate matrix G = F^-T -1 = (1+du)^-T - 1 = -[ (1+du)^-1 * du ]^T for each atom
     DO i=1, im
        Amat(:,:) = du(1:3,1:3,i)
        Amat(1,1) = Amat(1,1) + 1.d0
        Amat(2,2) = Amat(2,2) + 1.d0
        Amat(3,3) = Amat(3,3) + 1.d0
        CALL Mat3Inv(Amat, Bmat)
        !g(1:3,1:3,i) = Transpose(  Bmat(:,:) ) ! real Hartley & Mishin definition
        g(1:3,1:3,i) = -Transpose( MatMul( Bmat(:,:), du(1:3,1:3,i) ) ) ! Hartley & Mishin definition where we substract the identity matrix
     END DO

     ! Maximal number of neighbours
     nMax = MaxVal(nNeigh(1:im))

     ! Matrix containing neighbours positions
     ! ... in reference structure
     ALLOCATE(R0(1:3, 1:nMax))
     ! Matrix containing displacement gradient
     ALLOCATE(dG(1:nMax, 1:3, 1:3))

     ! Loop on atoms
     DO i=1, im
       nMax = nNeigh(i)
       ! Loop on neighbours
       DO n=1, nMax
          ! Index of neighbour atom
          j = iNeigh(n, i)
          R0(1:3,n) = xp0(1:3,j) - xp0(1:3,i)
          dG(n,1:3,1:3) = G(1:3,1:3,j) - G(1:3,1:3,i)
          IF (perioNyeTensor) THEN
                  ds0(:) = MatMul( inv_at0(:,:), R0(1:3,n) )
                  ds0(:) = aNInt(ds0(:))
                  R0(1:3,n) = R0(1:3,n) - MatMul( at0(:,:), ds0(:) )
          END IF
       END DO
       Bmat(:,:) = MatMul( R0(1:3,1:nMax), Transpose( R0(1:3,1:nMax) ) )
       CALL Mat3Inv(Bmat, invBmat)
       DO l=1, 3
          Vmat(:,:,l) = MatMul( R0(1:3,1:nMax), dG(1:nMax,1:3,l) )
          ! dGdx(i,j,l) = d G_j,l / d x_i
          dGdx(:,:,l) = MatMul( invBmat(:,:), Vmat(:,:,l) )
       END DO

       DO l=1, 3
          !!$ ! Initial version taken from Hartley and Mishin article
          !!$ ! nyeAtomic(1,l,i) = - d G_2,l / d x_3 + d G_3,l / d x_2 for atom i
          !!$ nyeAtomic(1,l,i) = dGdx(2,3,l) - dGdx(3,2,l)
          !!$ nyeAtomic(2,l,i) = dGdx(3,1,l) - dGdx(1,3,l)
          !!$ nyeAtomic(3,l,i) = dGdx(1,2,l) - dGdx(2,1,l)
          ! Correction of index error, as suggested by David Olmsted
          ! nyeAtomic(l,1,i) = + d G_2,l / d x_3 - d G_3,l / d x_2 for atom i
          nyeAtomic(l,1,i) = - dGdx(2,3,l) + dGdx(3,2,l)
          nyeAtomic(l,2,i) = - dGdx(3,1,l) + dGdx(1,3,l)
          nyeAtomic(l,3,i) = - dGdx(1,2,l) + dGdx(2,1,l)
       END DO
     END DO
     DEALLOCATE(R0) ; DEALLOCATE(dG)

  END SUBROUTINE BuildNyeTensor

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

  END SUBROUTINE mat3inv

END MODULE NyeTensorModule
