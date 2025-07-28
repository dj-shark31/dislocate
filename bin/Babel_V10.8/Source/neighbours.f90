MODULE NeighboursModule

  SAVE

  ! Distance between first neareast neighbours
  REAL(kind(0.d0)) :: rNeigh

  ! Use periodic boundary conditions to build neighbours table
  LOGICAL :: perioNeigh

  ! Maximal number of neighbours for each atom
  INTEGER :: max_nNeigh

  INTEGER, parameter, private :: verbosity_max=3

  PRIVATE :: mat3inv

CONTAINS

  SUBROUTINE BuildCell(at)

     IMPLICIT NONE
     ! Periodicity vectors at(1:3,1), at(1:3,2), ...
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at

     ! Lengths in each direction and corresponding direction
     REAL(kind(0.d0)) :: l1, l2, l3
     REAL(kind(0.d0)), dimension(1:3) :: n1Vec, n2Vec, n3Vec

     l1 = Sqrt( Sum( at(1:3,1)**2 ) )
     n1Vec(:) = at(:,1)/l1

     n2Vec(:) = at(:,2) - Sum( at(:,2)*n1Vec(:) )*n1Vec
     l2 = Sqrt( Sum( n2Vec(:)**2 ) )
     n2Vec(:) = n2Vec(:)/l2

     n3Vec(:) = at(:,3) - Sum( at(:,3)*n1Vec(:) )*n1Vec
     n3Vec(:) = n3Vec(:) - Sum( n3Vec(:)*n2Vec(:) )*n2Vec
     l3 = Sqrt( Sum( n3Vec(:)**2 ) )
     n3Vec(:) = n3Vec(:)/l3

  END SUBROUTINE BuildCell

  SUBROUTINE InitNeighbours(imm)

     IMPLICIT NONE

     ! Maximal number of atoms
     INTEGER, intent(in) :: imm

     perioNeigh=.true.

     IF (max_nNeigh.LE.0) max_nNeigh=16

  END SUBROUTINE InitNeighbours

  SUBROUTINE BuildNeighbours(im, xp, at, nNeigh, iNeigh, out)

     USE Babel_data, ONLY : verbosity
     IMPLICIT NONE

     ! Number of atoms
     INTEGER, intent(in) :: im
     ! Atoms coordinates
     REAL(kind(0.d0)), dimension(:,:), intent(in) :: xp
     ! Periodicity vectors at(1:3,1), at(1:3,2), ...
     !REAL(kind(0.d0)), dimension(1:3,1:3), intent(in), optional :: at
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
     ! Number of neighbours for each atom and maximal value
     INTEGER, dimension(:), intent(out) :: nNeigh
     ! Neighbour indexes for each atom
     INTEGER, dimension(:,:), intent(out) :: iNeigh

     ! Output unit 
     INTEGER, intent(in), optional :: out

     INTEGER :: i, j
     REAL(kind(0.d0)), dimension(1:3) :: dR, dS
     REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
     REAL(kind(0.d0)) :: rNeigh2, dR2


     IF (rNeigh.LE.0.d0) THEN
             WRITE(0,'(a)') "Distance between neighbours rNeigh needs to be positive"
             WRITE(0,'(a,g14.6)') "  rNeigh = ", rNeigh
             WRITE(0,'(a)') "Define rNeigh in input file"
             STOP "< BuildNeighbours >"
     END IF

     rNeigh2 = rNeigh**2

     IF (perioNeigh) THEN
             !IF (.NOT.Present(at)) THEN
                     !WRITE(0,'(a)') "Peridocitiy vectors need to be defined"
                     !STOP "< BuildNeighbours >"
             !END IF
             CALL Mat3Inv(at, inv_at)
     ELSE
             inv_at(:,:)=0.d0
     END IF

     DO i=1, im
        DO j=i+1, im
          ! Calculate distance between atoms i and j
          dR(:) = xp(:,j) - xp(:,i)
           IF (perioNeigh) THEN
                   ds(:) = MatMul( inv_at(:,:), dR(:) )
                   ds(:) = ds(:) - aNInt(ds(:))
                   dR(:) = MatMul( at(:,:), ds(:) )
           END IF
           dR2 = Sum( dR(:)**2 )

           IF (dR2.GT.rNeigh2) Cycle

           ! Atoms i and j are neighbours
           nNeigh(i) = nNeigh(i) + 1
           IF (nNeigh(i).GT.max_nNeigh) THEN
                   WRITE(0,'(a,i0)')    'Too many neighbours found for atom ', i
                   WRITE(0,'(a)')       'You need to increase max_nNeigh or to decrease rNeigh'
                   WRITE(0,'(a,i0)')    '  max_nNeigh = ', max_nNeigh
                   WRITE(0,'(a,g13.6)') '  rNeigh = ', rNeigh
                   STOP '< BuildNeighbours >'
           END IF
           iNeigh(nNeigh(i),i) = j

           nNeigh(j) = nNeigh(j) + 1
           IF (nNeigh(j).GT.max_nNeigh) THEN
                   WRITE(0,'(a,i0)')    'Too many neighbours found for atom ', j
                   WRITE(0,'(a)')       'You need to increase max_nNeigh or to decrease rNeigh'
                   WRITE(0,'(a,i0)')    '  max_nNeigh = ', max_nNeigh
                   WRITE(0,'(a,g13.6)') '  rNeigh = ', rNeigh
                   STOP '< BuildNeighbours >'
           END IF
           iNeigh(nNeigh(j),j) = i

        END DO
     END DO

     IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
             CALL PrintNeighbours(im, nNeigh, iNeigh, out)
     END IF

  END SUBROUTINE BuildNeighbours

  SUBROUTINE PrintNeighbours(im, nNeigh, iNeigh, out)

      IMPLICIT NONE

      ! Number of atoms
      INTEGER, intent(in) :: im
      ! Number of neighbours for each atom and maximal value
      INTEGER, dimension(:), intent(in) :: nNeigh
      ! Neighbour indexes for each atom
      INTEGER, dimension(:,:), intent(in) :: iNeigh
      ! Output unit
      INTEGER, intent(in) :: out

      INTEGER :: nMin, nMax, n, i

      nMin = MinVal(nNeigh(1:im)) ; nMax = MaxVal(nNeigh(1:im))

      WRITE(out,*)
      WRITE(out,'(a,g13.6)') 'Neighbour tables built for cutoff radius, rNeigh = ', rNeigh
      DO n=nMin, nMax
         i=Count( nNeigh(1:im).EQ.n )
         IF (i.GT.0) &
                 WRITE(out,'(2(a,i0))') '  number of atoms with ', n,' neigbours: ', i
      END DO
      WRITE(out,*)

  END SUBROUTINE PrintNeighbours

    SUBROUTINE mat3inv(A, B)
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

END MODULE NeighboursModule
