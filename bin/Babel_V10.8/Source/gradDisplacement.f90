MODULE GradDisplacementModule

   IMPLICIT NONE
        
   ! Use periodic boundary conditions to calculate displacement gradient ...
   LOGICAL, private :: perioStrainFromDisplacement

   ! Gradient of total and elastic displacement for each atom
   !  gradDisplacement(i,j, n) = d U_i / d X_j for atom n
   REAL(kind(0.d0)), dimension(:,:,:), allocatable :: gradDisplacement

  INTEGER, parameter, private :: verbosity_max=3
   ! Verbosity needed for printing debug information for calculation of the
   !   gradient of the total displacement
   INTEGER, parameter, private :: verbosity_debug_gradDisplacement=20

  PRIVATE :: mat3inv

CONTAINS

  SUBROUTINE InitGradDisplacement(imm)

     IMPLICIT NONE

     ! Maximal number of atoms
     INTEGER, intent(in) :: imm

     perioStrainFromDisplacement=.true.

     IF (Allocated(gradDisplacement)) Deallocate(gradDisplacement)
     ALLOCATE(gradDisplacement(1:3,1:3,1:imm))
     gradDisplacement(:,:,:) = 0.d0

  END SUBROUTINE InitGradDisplacement

  SUBROUTINE DestroyGradDisplacement

     IMPLICIT NONE

     IF (Allocated(gradDisplacement)) Deallocate(gradDisplacement)

  END SUBROUTINE DestroyGradDisplacement

  SUBROUTINE BuildGradDisplacement(im, xp0, u, at0, at, nNeigh, iNeigh, out)
     ! Calculate displacement gradient from finite differences using
     ! displacement of neighbours for each atom

     USE Babel_data, ONLY : verbosity
     IMPLICIT NONE

     ! Number of atoms
     INTEGER, intent(in) :: im
     ! Atoms coordinates in reference configuration
     REAL(kind(0.d0)), dimension(:,:), intent(in) :: xp0
     ! Displacement of each atom
     REAL(kind(0.d0)), dimension(:,:), intent(in) :: u
     ! Periodicity vectors of the reference structure: at0(1:3,1), at0(1:3,2), ...
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at0
     ! Periodicity vectors of the strained structure: at(1:3,1), at(1:3,2), ...
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
     ! Number of neighbours for each atom 
     INTEGER, dimension(:), intent(in) :: nNeigh
     ! Neighbour indexes for each atom i: iNeigh(1:nNeigh(i), i)
     INTEGER, dimension(:,:), intent(in) :: iNeigh
     ! Output unit 
     INTEGER, intent(in), optional :: out

     INTEGER :: i, j, nMax, n
     REAL(kind(0.d0)), dimension(1:3) :: ds0
     REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at0
     REAL(kind(0.d0)), dimension(:,:), allocatable :: R0, tR
     REAL(kind(0.d0)), dimension(1:3,1:3) :: Amat, Bmat, invBmat, Ftotal

     ! Inverse periodicity vectors
     IF (perioStrainFromDisplacement) THEN
             CALL Mat3Inv(at0, inv_at0)
     ELSE
             inv_at0(:,:)=0.d0
     END IF

     ! Maximal number of neighbours
     nMax = MaxVal(nNeigh(1:im))

     ! Matrix containing neighbours positions
     ! ... in reference structure
     ALLOCATE(R0(1:3, 1:nMax))
     ! ... in strained structure (transposed matrix)
     ALLOCATE(tR(1:nMax, 1:3))

     ! Loop on atoms
     DO i=1, im
       nMax = nNeigh(i)
       ! Loop on neighbours
       DO n=1, nMax
          ! Index of neighbour atom
          j = iNeigh(n, i)
          R0(1:3,n) = xp0(1:3,j) - xp0(1:3,i)
          tR(n,1:3) = u(1:3,j) - u(1:3,i)
          IF (perioStrainFromDisplacement) THEN
                  ds0(:) = MatMul( inv_at0(:,:), R0(1:3,n) )
                  ds0(:) = aNInt(ds0(:))
                  R0(1:3,n) = R0(1:3,n) - MatMul( at0(:,:), ds0(:) )
                  tR(n,1:3) = tR(n,1:3) - MatMul( at(:,:)-at0(:,:), ds0(:) )
          END IF
       END DO
       Amat(:,:) = MatMul( R0(1:3,1:nMax), tR(1:nMax,1:3) )
       Bmat(:,:) = MatMul( R0(1:3,1:nMax), Transpose( R0(1:3,1:nMax) ) )
       CALL Mat3Inv(Bmat, invBmat)
       ! Ftotal(i,j) = d U_j / d X_i
       Ftotal(:,:) = MatMul( Amat(:,:), invBmat(:,:) )
       ! gradDisplacement(i,j, n) = d U_i / d X_j for atom n
       gradDisplacement(1:3,1:3,i) = Transpose( Ftotal(:,:) )

       IF (verbosity.GE.verbosity_debug_gradDisplacement) THEN
               WRITE(out,'(a,i0)') ' atom ', i
               DO n=1, nMax
                  WRITE(out,'(a,i3,2(a,3g14.6))') '  neighbour ', n, &
                          ' R(1:3) = ',  R0(1:3,n)+tR(n,1:3), &
                          ' - R0(1:3) = ', R0(1:3,n)
               END DO
               WRITE(out,'(a,3g14.6)') '  dUx / dR(1:3) = ', Ftotal(1,1:3)
               WRITE(out,'(a,3g14.6)') '  dUy / dR(1:3) = ', Ftotal(2,1:3)
               WRITE(out,'(a,3g14.6)') '  dUz / dR(1:3) = ', Ftotal(3,1:3)
               WRITE(out,*)
       END IF

     END DO

     DEALLOCATE(R0) ; DEALLOCATE(tR)
        
  END SUBROUTINE BuildGradDisplacement

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

END MODULE GradDisplacementModule
