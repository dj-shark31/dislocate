PROGRAM Burgers2D

  !USE Math
  IMPLICIT NONE

  REAL(kind(0.d0)), parameter :: Zero=1.d-6
  REAL(kind(0.d0)), parameter :: distance_zero=1.d-1

  CHARACTER(len=200) :: inpFile
  INTEGER :: i, im, n, i0Up, i1Up, i0Down, i1Down, nUp, nDown, &
         nUpLines, nDownLines, n0UpLine, n0DownLine
  INTEGER, dimension(:), allocatable :: iUp, iDown, UpLines, DownLines
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp, u  
  REAL(kind(0.d0)), dimension(:,:), allocatable :: sUp, sDown
  REAL(kind(0.d0)), dimension(:), allocatable :: yUpLine, yDownLine
  REAL(kind(0.d0)), dimension(1:3) :: x0, n0, lDislo, lScan, uPerio
  REAL(kind(0.d0)) :: d, d0, half_d0
  CHARACTER(len=5) :: label_temp
  CHARACTER(len=11) :: outFile
  LOGICAL :: first


  ! Read name of the input file containing atom coordinates and displacements
  WRITE(6,'(a)') '# Name of the input file with atom displacements (xyz format)'
  READ(5,*) inpFile

  ! Read input structure file
  OPEN(file=inpFile, unit=51, action="read", status="old")
  READ(51,*) im         ! Number of atoms
  READ(51,*)            ! Line for title
  ALLOCATE(xp(1:3,1:im)) 
  ALLOCATE(u(1:3,1:im))
  DO i=1, im
     READ(51,*) label_temp, xp(1:3,i), u(1:3,i)
  END DO
  CLOSE(51)

  ! Read definition of the cut plane
  WRITE(6,'(a)') '# Coordinates of a point belonging to the cut plane'
  READ(5,*) x0(1:3)
  WRITE(6,'(a)') '# Normal to the cut plane'
  READ(5,*) n0(1:3)
  n0(1:3) = n0(1:3)/Sqrt( Sum( n0(1:3)**2 ) )
  WRITE(6,'(a)') '# Interplane distance'
  READ(5,*) d0
  half_d0=0.51*d0

  ! Definition of the dislocation
  WRITE(6,'(a)') '# Line direction of the dislocation'
  READ(5,*) lDislo(1:3)
  lDislo(1:3) = lDislo(1:3)/Sqrt( Sum( lDislo(1:3)**2 ) )

  ! Check that dislocation line direction is perpendicular to cut plane normal
  d = Sum( lDislo(1:3)*n0(1:3) )
  IF (Abs(d).GT.Zero) THEN
          WRITE(0,'(a)') 'Dislocation line direction is not perpendicular to cut plane normal'
          STOP
  END IF

  ! Third direction is perpendicular to dislocation line and to cut plane normal
  lScan(1:3) = CrossProduct( lDislo(1:3), n0(1:3) )

  ! Use of periodic conditions to calculate displacement
  WRITE(6,'(a)') '# Shortest periodicity length in each direction' 
  WRITE(6,'(a)') '#  (<= 0 for not applying periodic conditions on displacement)'
  READ(5,*) uPerio(1:3)

  ! Look for atoms just above and below the cut plane
  ALLOCATE(iUp(1:im)) ; ALLOCATE(iDown(1:im))
  iUp(:)=0 ; iDown(:)=0
  nUp=0 ; nDown=0
  DO i=1, im
     ! Distance of the atom to the glide plane
     d = Sum( ( xp(1:3,i)-x0(1:3) )*n0(1:3) )
     IF ( (d.GE.0.d0).AND.(d.LE.half_d0) ) THEN
             ! Atom belongs to upper plane
             nUp = nUp + 1
             iUp(nUp) = i
     ELSE IF ( (d.LE.0.d0).AND.(d.GE.-half_d0) ) THEN
             ! Atom belong to lower plane
             nDown = nDown + 1
             iDown(nDown) = i
     END IF
  END DO


  ! Coordinates along the direction iScan(1:3)
  IF (nUp.EQ.0) THEN
          WRITE(0,'(a)') 'Did not find any atom in the upper plane'
          WRITE(0,'(a,g14.6)') 'Check value of interplane distance: d0 = ', d0
          STOP '< mainBurgersDistri >'
  END IF
  ALLOCATE(sUp(1:2,1:nUp))
  DO n=1, nUp
     sUp(1,n) = Sum( ( xp(1:3,iUp(n)) - x0(1:3) )*lScan(1:3) )
     sUp(2,n) = Sum( ( xp(1:3,iUp(n)) - x0(1:3) )*lDislo(1:3) )
  END DO

  IF (nDown.EQ.0) THEN
          WRITE(0,'(a)') 'Did not find any atom in the lower plane'
          WRITE(0,'(a,g14.6)') 'Check value of interplane distance: d0 = ', d0
          STOP '< mainBurgersDistri >'
  END IF
  ALLOCATE(sDown(1:2,1:nDown))
  DO n=1, nDown
     sDown(1,n) = Sum( ( xp(1:3,iDown(n)) - x0(1:3) )*lScan(1:3) )
     sDown(2,n) = Sum( ( xp(1:3,iDown(n)) - x0(1:3) )*lDislo(1:3) )
  END DO

  ! Sort atoms according to y
  CALL TriAtoms(sUp, iUp, nUp, 2)
  CALL TriAtoms(sDown, iDown, nDown, 2)

  ! Gather atoms in lines sharing the ~ same y
  ALLOCATE(UpLines(1:nUp))
  ALLOCATE(DownLines(1:nDown))
  CALL GatherLines(sUp, nUp, 2, nUpLines, UpLines)
  CALL GatherLines(sDown, nDown, 2, nDownLines, DownLines)

  ! Sort atoms according to x in each line 
  ! and create temporary table containing y coordinate of each line
  Allocate(yUpLine(1:nUpLines))
  i0Up=0
  DO n=1, nUpLines
     i1Up=i0Up+UpLines(n)
     CALL TriAtoms(sUp(:,i0Up+1:i1Up), iUp(i0Up+1:i1Up), i1Up-i0Up, 1)
     yUpLine(n)=sUp(2,i0Up+1)
     i0Up=i1Up
  END DO

  Allocate(yDownLine(1:nDownLines))
  i0Down=0
  DO n=1, nDownLines
     i1Down=i0Down+DownLines(n)
     CALL TriAtoms(sDown(:,i0Down+1:i1Down), iDown(i0Down+1:i1Down), i1Down-i0Down, 1)
     yDownLine(n)=sDown(2,i0Down+1)
     i0Down=i1Down
  END DO

  ! Look for the closest starting lines on upper and lower planes
  CALL FindStartingPoint(yUpLine, yDownLine, n0UpLine, n0DownLine, first=.true.)
 
  WRITE(6,*)
  WRITE(6,'(2a)') '# Name of the input file with atom displacements: ', &
        Trim(inpFile)
  WRITE(6,'(a,3g14.6)') '# Coordinates of a point belonging to the cut plane: x0 = ', &
        x0(1:3)
  WRITE(6,'(a,3g14.6)') '# Normal to the cut plane: n = ', n0(1:3)
  WRITE(6,'(a,g14.6)') '# Interplane distance: d = ', d0
  WRITE(6,'(a,3g14.6)') '# Line direction of the dislocation: l = ', lDislo(1:3)
  WRITE(6,'(a,3g14.6)') '# Shortest periodicity length in each direction: ', &
        uPerio(1:3)
  WRITE(6,*)
  WRITE(6,'(a,3g14.6)') '# s(1): point projected in the glide plane along l^n = ', lScan(1:3)
  WRITE(6,'(a,3g14.6)') '# s(2): point projected in the glide plane along  l  = ', lDislo(1:3)
  WRITE(6,'(a,3g14.6)') '# u(1): displacement difference projected along x = ', &
        1.d0, 0.d0, 0.d0
  WRITE(6,'(a,3g14.6)') '# u(2): displacement difference projected along y = ', &
        0.d0, 1.d0, 0.d0
  WRITE(6,'(a,3g14.6)') '# u(3): displacement difference projected along z = ', &
        0.d0, 0.d0, 1.d0
  WRITE(6,*)

  i0Down=0
  DO n=1, n0DownLine      ! Loop on y-lines
     ! Consider line n for Down plane
     i1Down=i0Down+DownLines(n)
     i0Down=i1Down
  END DO

  i0Up=0
  DO n=1, n0UpLine      ! Loop on y-lines
     ! Consider line n for upper plane
     i1Up=i0Up+UpLines(n)
     i0Up=i1Up
  END DO

  DO n=1, Min(nUpLines-n0UpLine, nDownLines-n0DownLine)
     ! Consider lines: n+n0UpLine for upper plane
     !                 n+n0DownLine for down plane
     i1Up=i0Up+UpLines(n+n0UpLine)
     i1Down=i0Down+DownLines(n+n0DownLine)
     IF (n.EQ.1) THEN
             first=.true.
     ELSE
             first=.FALSE.
     END IF

     ! Write on the screen
     WRITE(6,'(3(a,i0),a)') "# Line ", n, " containing ", &
        i1Up-i0Up, " points in upper plane and ", &
        i1Down-i0Down, " points in lower plane"
     CALL PrintLine( sUp(:,i0Up+1:i1Up), iUp(i0Up+1:i1Up), &
        sDown(:,i0Down+1:i1Down), iDown(i0Down+1:i1Down), 6, first )

     ! Write in a file "line???.res"
     IF (n.LE.9) THEN
             WRITE(outFile,'(a6,i1,a4)') "line00", n, ".res"
     ELSE IF (n.LE.99) THEN
             WRITE(outFile,'(a5,i2,a4)') "line0", n, ".res"
     ELSE IF (n.LE.999) THEN
             WRITE(outFile,'(a4,i3,a4)') "line", n, ".res"
     ELSE
             WRITE(0,'(a)') "Line index n cannot be bigger than 999"
             STOP
     END IF
     OPEN(file=outFile,unit=60,action='write')
     WRITE(60,'(3(a,i0),a)') "# Line ", n, " containing ", &
        i1Up-i0Up, " points in upper plane and ", &
        i1Down-i0Down, " points in lower plane"
     CALL PrintLine( sUp(:,i0Up+1:i1Up), iUp(i0Up+1:i1Up), &
        sDown(:,i0Down+1:i1Down), iDown(i0Down+1:i1Down), 60, first )
      CLOSE(60)

     i0Up=i1Up
     i0Down=i1Down
     WRITE(6,*)
  END DO


  DEALLOCATE(yUpLine) ; DEALLOCATE(yDownLine)
  DEALLOCATE(UpLines) ; DEALLOCATE(DownLines)
  DEALLOCATE(sUp) ; DEALLOCATE(sDown)
  DEALLOCATE(iUp) ; DEALLOCATE(iDown)
  DEALLOCATE(xp) ; DEALLOCATE(u)

CONTAINS 

  !================================================
  
  SUBROUTINE TriAtoms(s, iAtom, nPoints, i0)
    ! Sort positions of atoms according to coordinate i0

    IMPLICIT NONE
    ! Atom coordinates
    REAL(kind(0.d0)), dimension(:,:), intent(inout) :: s
    ! Atom index
    INTEGER, dimension(:), intent(inout) :: iAtom
    ! Number of atoms
    INTEGER, intent(in) :: nPoints
    ! Coordinate index to sort
    INTEGER, intent(in) :: i0

    INTEGER :: n, m, i_temp
    REAL(kind(0.d0)), dimension(1:Size(s,1)) :: s_temp

    DO n=1, nPoints
       m = MinLoc(s(i0,n:nPoints), 1) + n-1
       s_temp(:) = s(:,n) ; i_temp = iAtom(n)
       s(:,n) = s(:,m) ; iAtom(n) = iAtom(m)
       s(:,m) = s_temp(:) ; iAtom(m) = i_temp
    END DO

  END SUBROUTINE TriAtoms

  !================================================

  SUBROUTINE GatherLines(s, nPoints, i0, nLines, lines)
    ! Gather atoms in lines sharing same coordinate i0

    IMPLICIT NONE
    ! Atom coordinates
    REAL(kind(0.d0)), dimension(:,:), intent(in) :: s
    ! Number of atoms
    INTEGER, intent(in) :: nPoints
    ! Coordinate index to sort
    INTEGER, intent(in) :: i0
    ! Number of different lines
    INTEGER, intent(out) :: nLines
    ! Number of atoms in each lines
    INTEGER, dimension(:), intent(out) :: lines

    REAL(kind(0.d0)) :: sOld
    INTEGER :: n

    nLines=0
    Lines(:)=0

    IF (nPoints.LE.0) RETURN

    ! Initialization for first point
    nLines=1
    lines(1)=1
    sOld=s(i0,1)

    DO n=2, nPoints
       IF ( Abs(s(i0,n)-sOld).LE.distance_zero ) THEN
               ! Same line
               lines(nLines)=lines(nLines)+1
       ELSE
               ! New line
               nLines=nLines+1
               lines(nLines)=1
       END IF
       sOld=s(i0,n)
    END DO

  END SUBROUTINE GatherLines

  !================================================

  SUBROUTINE FindStartingPoint(sUp, sDown, n0Up, n0Down, first)
  ! Look for the closest starting points on upper and lower planes

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:), intent(in) :: sUp, sDown
    INTEGER, intent(out) :: n0Up, n0Down
    LOGICAL, intent(in) :: first

    REAL(kind(0.d0)) :: s0Up, s0Down
    INTEGER :: nUp, nDown
    REAL(kind(0.d0)), save :: s0Up_moins_s0Down

    nUp=Size(sUp, 1) ; nDown=Size(sDown, 1)

    s0Up = sUp(1) ; s0Down = sDown(1)
    n0Up=1 ; n0Down=1

    IF ( (s0Up.GT.s0Down).AND.(nDown.GT.1) ) THEN
            n0Up = 1
            n0Down = 2
            DO WHILE (sDown(n0Down).LT.s0Up)
               IF (n0Down+1.LE.nDown) THEN
                       n0Down = n0Down+1
               ELSE
                       EXIT
               END IF
            END DO
            ! We should have sDown(n0Down-1) < sUp(n0Up) < sDown(n0Down)
            IF (first) THEN
                    IF ((s0Up-sDown(n0Down-1)).LT.(sDown(n0Down)-s0Up)) n0Down=n0Down-1
                    s0Up_moins_s0Down = sUp(n0Up) - sDown(n0Down)
            ELSE
                    IF (s0Up_moins_s0Down.GT.0.d0) n0Down=n0Down-1
            END IF
    ELSE IF ( (s0Down.GT.s0Up).AND.(nUp.GT.1) ) THEN
            n0Down = 1
            n0Up = 2
            DO WHILE (sUp(n0Up).LT.s0Down)
               IF (n0Up+1.LE.nUp) THEN
                       n0Up = n0Up+1
               ELSE
                       EXIT
               END IF
            END DO
            ! We should have sUp(n0Up-1) < sDown(n0Down) < sUp(n0Up)
            IF (first) THEN
                    IF ((s0Down-sUp(n0Up-1)).LT.(sUp(n0Up)-s0Down)) n0Up=n0Up-1
                    s0Up_moins_s0Down = sUp(n0Up) - sDown(n0Down)
            ELSE
                    IF ( s0Up_moins_s0Down.LT.0.d0 ) n0Up=n0Up-1
            END IF
    END IF
    n0Up=n0Up-1 ; n0Down=n0Down-1

  END SUBROUTINE FindStartingPoint

  !================================================

  SUBROUTINE PrintLine(sUp, iUp, sDown, iDown, out, first)

    IMPLICIT NONE

    REAL(kind(0.d0)), dimension(:,:), intent(in) :: sUp, sDown
    INTEGER, dimension(:), intent(in) :: iUp, iDown
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: first

    INTEGER :: n0Up, n0Down, n, i, nUp, nDown
    REAL(kind(0.d0)), dimension(1:3) :: du
    REAL(kind(0.d0)), dimension(1:2) :: s

    nUp=Size(sUp, 2) ; nDown=Size(sDown, 2)

    ! Look for the closest starting points on upper and lower planes
    CALL FindStartingPoint(sUp(1,:), sDown(1,:), n0Up, n0Down, first)

    ! Write results on output unit
    WRITE(out,'(a)') '# s(1:2), u(1:3), xDown(1:3), xUp(1:3)'
    IF (n0Down.GT.0) THEN
            DO n=1, n0Down
               WRITE(out,'(5a14,3g14.6,3a14)') &
                  '?', '?', '?', '?', '?', xp(1:3,iDown(n)), '?', '?', '?'
            END DO
    END IF
    IF (n0UP.GT.1) THEN
            DO n=1, n0Up
               WRITE(out,'(5a14,3a14,3g14.6)') &
                  '?', '?', '?', '?', '?', '?', '?', '?', xp(1:3,iUp(n))
            END DO
    END IF

    !!$duOld(:)=0.d0
    DO n=1, Min(nUp-n0Up, nDown-n0Down)
       s(1:2) = 0.5d0*( sUp(1:2,n+n0Up) + sDown(1:2,n+n0Down) )
       du(1:3) = u(1:3,iUp(n+n0Up)) - u(1:3,iDown(n+n0Down))               
       ! Apply periodic boundary conditions 
       DO i=1, 3
          IF (uPerio(i).LE.distance_zero) Cycle
          du(i) = Modulo( du(i), uPerio(i) )
       END DO
       WRITE(out,'(11g14.6)') s(1:2), du(1:3), xp(1:3,iDown(n+n0Down)), xp(1:3,iUp(n+n0Up))
    END DO

    IF (nDown-n0Down.GT.nUp-n0Up) THEN
            DO n=nUp-n0Up+n0Down+1, nDown
               WRITE(out,'(5a14,3g14.6,3a14)') &
                  '?','?','?','?','?',xp(1:3,iDown(n)), '?', '?', '?'
            END DO
    ELSE IF (nUP-n0up.GT.nDown-n0Down) THEN
            DO n=nDown-n0Down+n0Up+1, nUp
               WRITE(out,'(8a14,3g14.6)') &
                  '?','?','?','?','?','?', '?', '?', xp(1:3,iUp(n))
            END DO
    END IF


  END SUBROUTINE PrintLine

  !=========================================

  FUNCTION CrossProduct(ex,ey) RESULT(ez)
    
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: ex, ey
    REAL(kind(0.d0)), dimension(1:3) :: ez

    ez(1) = ex(2)*ey(3) - ex(3)*ey(2)
    ez(2) = ex(3)*ey(1) - ex(1)*ey(3)
    ez(3) = ex(1)*ey(2) - ex(2)*ey(1)

  END FUNCTION CrossProduct

END PROGRAM Burgers2D
