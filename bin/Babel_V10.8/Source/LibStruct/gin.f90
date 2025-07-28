MODULE gin_module
  ! Structure files in gin format (compatible with NDM90)

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter, private :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter, private :: Distance_Zero2=Distance_Zero*Distance_Zero
 
CONTAINS

  !=============================================================
  
  FUNCTION GetImmGin(inp) RESULT(imm)
    ! Get maximal number of atoms imm from gin file connected to unit inp
    
    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    INTEGER :: imm

    INTEGER :: xLat, yLat, zLat, imcell

    imm=0

    ! Number of cells in 3 directions
    READ(inp,*) xLat, yLat, zLat
    IF ( (xLat.LT.1).OR.(yLat.LT.1).OR.(zLat.LT.1) ) THEN
            WRITE(0,'(a)') 'Problem with number of times the unit cell has to be duplicated'
            WRITE(0,'(a,3(i0,1x),a)') ' lat(1:3) = ', xLat, yLat, zLat, ' (read in gin file)'
            STOP '< GetImmGin >'
    END IF
        
    ! Lattice vector coordinates of the cell
    READ(inp,*) 
    READ(inp,*)
    READ(inp,*)

    ! Number of atoms in the cell
    READ(inp,*) imcell

    ! Number of atoms in the simulation box
    imm=imcell*xLat*yLat*zLat


    REWIND(inp)

  END FUNCTION GetImmGin

  !=============================================================
  
  SUBROUTINE ReadGin(xp, iTyp, im, at, nTypes, inp)

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(out) ::xp
    INTEGER, dimension(:), intent(out) :: iTyp
    INTEGER, intent(out) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: at
    INTEGER, intent(inout) :: nTypes
    INTEGER, intent(in) :: inp

    INTEGER :: i, j,k,l, j0, n
    ! Number of atoms in unit cell
    INTEGER :: imCell, xLat, yLat, zLat
    REAL(kind(0.d0)), dimension(:,:), allocatable :: xc

    ! Number of cells in 3 directions
    READ(inp,*) xLat, yLat, zLat
    IF ( (xLat.LT.1).OR.(yLat.LT.1).OR.(zLat.LT.1) ) THEN
            WRITE(0,'(a)') 'Problem with number of times the unit cell have to be duplicated'
            WRITE(0,'(a,3(i0,1x),a)') ' lat(1:3) = ', xLat, yLat, zLat, ' (read in gin file)'
            STOP '< ReadGin >'
    END IF
    
    ! Lattice vector coordinates of the cell
    READ(inp,*) at(1:3,1)
    READ(inp,*) at(1:3,2)
    READ(inp,*) at(1:3,3)

    ! Number of atoms in the cell
    READ(inp,*) imcell

    ! Number of atoms in the simulation box
    im=imcell*xLat*yLat*zLat

    ! Maximal number of atoms in simulation box (used to allocate tables)
    ! Maximal number of atoms in simulation box
    IF ( (im.GT.size(xp,2)).OR.(im.GT.size(iTyp,1)) ) THEN
             WRITE(0,'(a,i0)') '  dimension of xp(1:3,:): ', size(xp,2)
             WRITE(0,'(a,i0)') '  dimension of iTyp(:): ', size(iTyp,1)
             WRITE(0,'(a,i0)') '  number of atoms read in GIN file: ', im
             STOP '< ReadGin >'
    END IF

    IF (Allocated(xc)) Deallocate(xc)
    ALLOCATE(xc(1:3,1:imcell))

    DO i=1, imcell
      READ(inp,*) xc(1:3,i), ityp(i)
    END DO

    ! Atom real coordinates
    xp(1:3,1:imcell) = MatMul( at(1:3,1:3), xc(1:3,1:imcell) )
    Deallocate(xc)

    ! Duplicate unit cell
    DO l=1, zLat
       DO k=1, yLat
          IF (k*l.EQ.1) THEN
                  j0=2
          ELSE
                  j0=1
          END IF
          DO j=j0, xLat
             DO i=1, imcell
                n = i + imcell*( (l-1)*yLat*xLat + (k-1)*xLat +j-1 )
                xp(1:3,n ) = xp(1:3,i) &
                   + dble(j-1)*at(1:3,1) + dble(k-1)*at(1:3,2) + dble(l-1)*at(1:3,3)
                ityp(n) = ityp(i)
             END DO
          END DO
       END DO
    END DO

    ! Lattice vector coordinates of the simulation box
    at(1:3,1) = dble(xLat)*at(1:3,1)
    at(1:3,2) = dble(yLat)*at(1:3,2)
    at(1:3,3) = dble(zLat)*at(1:3,3)

    ! Number of atom types
    n=MaxVal(ityp(1:im))
    nTypes=Max(nTypes,n)

  END SUBROUTINE ReadGin
 
  !=============================================================

  SUBROUTINE WriteGin(xp, iTyp, im, at, out, mask)
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written

    USE MathStruct
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
    INTEGER, dimension(:), intent(in) :: iTyp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    INTEGER, intent(in) :: out
    LOGICAL, dimension(:), intent(in), optional :: mask

    ! Local variables
    INTEGER :: i
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
    REAL(kind(0.d0)), dimension(:,:), allocatable :: xc

    IF ( SUM( at(1:3,1:3)**2 ).LE.Distance_Zero2 ) THEN
            WRITE(0,'(a)') 'You need to define lattice vectors at(1:3,i) to be&
                & able to use gin format'
            STOP '< WriteGin >'
    END IF
    CALL Mat3Inv(at,inv_at)
    IF (Allocated(xc)) Deallocate(xc)
    ALLOCATE(xc(1:3,1:im))
    xc(1:3,1:im) = MatMul( inv_at(1:3,1:3), xp(1:3,1:im) )

    WRITE(out,'(3(i0,1x))') 1, 1, 1
    WRITE(out,'(3(g24.16,1x))') at(1:3,1)
    WRITE(out,'(3(g24.16,1x))') at(1:3,2)
    WRITE(out,'(3(g24.16,1x))') at(1:3,3)
    IF (Present(mask)) THEN
            WRITE(out,'(i0)') Count(mask(1:im))
            DO i=1, im
               IF (mask(i)) WRITE(out,'(3g24.16,2x,i0)') xc(1:3,i), ityp(i)
            END DO
    ELSE
            WRITE(out,'(i0)') im
            DO i=1, im
              WRITE(out,'(3g24.16,2x,i0)') xc(1:3,i), ityp(i)
            END DO
    END IF

    Deallocate(xc)

  END SUBROUTINE WriteGin

END MODULE gin_module
