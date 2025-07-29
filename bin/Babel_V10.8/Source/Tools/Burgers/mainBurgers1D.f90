PROGRAM Burgers1D

  IMPLICIT NONE

  REAL(kind(0.d0)), parameter :: Zero=1.d-6
  REAL(kind(0.d0)), parameter :: distance_zero=1.d-4

  CHARACTER(len=200) :: inpFile
  INTEGER :: i, im, n, m, nUp, nDown, i_temp, n0Up, n0Down
  INTEGER, dimension(:), allocatable :: iUp, iDown
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp, u  
  REAL(kind(0.d0)), dimension(:), allocatable :: sUp, sDown
  REAL(kind(0.d0)), dimension(1:3) :: cCut, nCut, dx, lLine, lScan, du, uPerio, duOld
  REAL(kind(0.d0)), dimension(1:3) :: uUp, uDown
  REAL(kind(0.d0)) :: d, dCut, half_dCut, s_temp, s, s0Up, s0Down
  INTEGER :: dnUp
  LOGICAL :: new_atom
  CHARACTER(len=5) :: label_temp


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
  WRITE(6,'(a)') '# Coordinates of a point belongin to the cut plane'
  READ(5,*) cCut(1:3)
  WRITE(6,'(a)') '# Normal to the cut plane'
  READ(5,*) nCut(1:3)
  nCut(1:3) = nCut(1:3)/Sqrt( Sum( nCut(1:3)**2 ) )
  WRITE(6,'(a)') '# Interplane distance'
  READ(5,*) dCut
  half_dCut=0.51*dCut

  ! Definition of the dislocation
  WRITE(6,'(a)') '# Line direction of the dislocation'
  READ(5,*) lLine(1:3)
  lLine(1:3) = lLine(1:3)/Sqrt( Sum( lLine(1:3)**2 ) )

  ! Check that dislocation line direction is perpendicular to cut plane normal
  d = Sum( lLine(1:3)*nCut(1:3) )
  IF (Abs(d).GT.Zero) THEN
          WRITE(0,'(a)') 'Dislocation line direction is not perpendicular to cut plane normal'
          STOP
  END IF

  ! Third direction is perpendicular to dislocation line and to cut plane normal
  lScan(1:3) = CrossProduct( lLine(1:3), nCut(1:3) )

  ! Use of periodic conditions to calculate displacement
  WRITE(6,'(a)') '#Shortest periodicity length in each direction' 
  WRITE(6,'(a)') '#  (<= 0 for not applying periodic conditions on displacement)'
  READ(5,*) uPerio(1:3)

  ! Index offset for selecting atoms in upper plane
  WRITE(6,'(a)') "Index offset for selecting atoms in upper plane  (default: 0)"
  READ(5,*) dnUp

  ! Look for atoms just above and below the cut plane
  ALLOCATE(iUp(1:im)) ; ALLOCATE(iDown(1:im))
  iUp(:)=0 ; iDown(:)=0
  nUp=0 ; nDown=0
  DO i=1, im
     d = Sum( ( xp(1:3,i)-cCut(1:3) )*nCut(1:3) )
     IF ( (d.GE.0.d0).AND.(d.LE.half_dCut) ) THEN
             ! Atom belongs to upper plane
             ! Check if this atom is not a replica of another atom
             new_atom=.TRUE.
             DO n=1, nUp
                dx(1:3) = xp(1:3,i) - xp(1:3,iUp(n))
                d = Sum( dx(1:3)*lScan(1:3) )
                IF (Abs(d).LE.distance_zero) THEN
                        new_atom=.FALSE.
                        exit
                END IF
             END DO
             IF (new_atom) THEN
                     ! We have a new atom
                     nUp = nUp + 1
                     iUp(nUp) = i
             ELSE
                     ! Chosse which atom to keep
                     d = Sum( dx(1:3)*lLine(1:3) )
                     IF (d.LT.0.d0) iUp(n) = i
             END IF
     ELSE IF ( (d.LE.0.d0).AND.(d.GE.-half_dCut) ) THEN
             ! Atom belong to lower plane
             ! Check if this atom is not a replica of another atom
             new_atom=.TRUE.
             DO n=1, nDown
                dx(1:3) = xp(1:3,i) - xp(1:3,iDown(n))
                d = Sum( dx(1:3)*lScan(1:3) )
                IF (Abs(d).LE.distance_zero) THEN
                        new_atom=.FALSE.
                        exit
                END IF
             END DO
             IF (new_atom) THEN
                     ! We have a new atom
                     nDown = nDown + 1
                     iDown(nDown) = i
             ELSE
                     ! Choose which atom to keep
                     d = Sum( dx(1:3)*lLine(1:3) )
                     IF (d.LT.0.d0) iDown(n) = i
             END IF
     END IF
  END DO


  ! Coordinates along the direction iScan(1:3)
  IF (nUp.EQ.0) THEN
          WRITE(0,'(a)') 'Did not find any atom in the upper plane'
          WRITE(0,'(a,g18.6)') 'Check value of interplane distance: dCut = ', dCut
          STOP '< mainBurgersDistri >'
  END IF
  ALLOCATE(sUp(1:nUp))
  DO n=1, nUp
     sUp(n) = Sum( ( xp(1:3,iUp(n)) - cCut(1:3) )*lScan(1:3) )
  END DO

  IF (nDown.EQ.0) THEN
          WRITE(0,'(a)') 'Did not find any atom in the lower plane'
          WRITE(0,'(a,g18.6)') 'Check value of interplane distance: dCut = ', dCut
          STOP '< mainBurgersDistri >'
  END IF
  ALLOCATE(sDown(1:nDown))
  DO n=1, nDown
     sDown(n) = Sum( ( xp(1:3,iDown(n)) - cCut(1:3) )*lScan(1:3) )
  END DO

  ! Sort atoms
  DO n=1, nUp
     m = MinLoc(sUp(n:nUp), 1) + n-1
     s_temp = sUp(n) ; i_temp = iUp(n)
     sUp(n) = sUp(m) ; iUp(n) = iUp(m)
     sUp(m) = s_temp ; iUp(m) = i_temp
  END DO
  DO n=1, nDown
     m = MinLoc(sDown(n:nDown), 1) + n-1
     s_temp = sDown(n) ; i_temp = iDown(n)
     sDown(n) = sDown(m) ; iDown(n) = iDown(m)
     sDown(m) = s_temp ; iDown(m) = i_temp
  END DO

 
  ! Look for the closest starting points on upper and lower planes
  s0Up = sUp(1) ; s0Down = sDown(1)
  IF (s0Up.GT.s0Down) THEN
          n0Up = 1
          n0Down = 2
          DO WHILE (sDown(n0Down).LT.s0Up)
             n0Down = n0Down+1
          END DO
          ! We should have sDown(n0Down-1) < sUp(n0Up) < sDown(i0)
          IF ((s0Up-sDown(n0Down-1)).LT.(sDown(n0Down)-s0Up)) n0Down=n0Down-1
  ELSE
          n0Down = 1
          n0Up = 2
          DO WHILE (sUp(n0Up).LT.s0Down)
             n0Up = n0Up+1
          END DO
          ! We should have sUp(n0Up-1) < sDown(n0Down) < sUp(i0)
          IF ((s0Down-sUp(n0Up-1)).LT.(sUp(n0Up)-s0Down)) n0Up=n0Up-1
  END IF
  n0Up=n0Up-1 ; n0Down=n0Down-1 

  IF (dnUp.GT.0) THEN
          n0Up=n0Up+dnUp
  ELSE
          n0Down=n0Down-dnUp
  END IF

  WRITE(6,*)
  WRITE(6,'(2a)') '# Name of the input file with atom displacements: ', &
        Trim(inpFile)
  WRITE(6,'(a,3g14.6)') '# Coordinates of a point belonging to the cut plane: cCut = ', &
        cCut(1:3)
  WRITE(6,'(a,3g14.6)') '# Normal to the cut plane: n = ', nCut(1:3)
  WRITE(6,'(a,g14.6)') '# Interplane distance: d = ', dCut
  WRITE(6,'(a,3g14.6)') '# Line direction of the dislocation: l = ', lLine(1:3)
  WRITE(6,'(a,3g14.6)') '# Shortest periodicity length in each direction: ', &
        uPerio(1:3)
  WRITE(6,'(a,i0)') "Index offset for selecting atoms in upper plane: ", dnUp
  WRITE(6,*)
  WRITE(6,'(a,3g14.6)') '# s(1): point projected in the glide plane along l^n = ', lScan(1:3)
  WRITE(6,'(a,3g14.6)') '# u(1): displacement difference projected along x = ', &
        1.d0, 0.d0, 0.d0
  WRITE(6,'(a,3g14.6)') '# u(2): displacement difference projected along y = ', &
        0.d0, 1.d0, 0.d0
  WRITE(6,'(a,3g14.6)') '# u(3): displacement difference projected along z = ', &
        0.d0, 0.d0, 1.d0
  WRITE(6,*)


  ! Write results on output unit
  WRITE(6,*)
  WRITE(6,'(a)') '# s, u(1:3), xDown(1:3), xUp(1:3), sDown, sUp'
  IF (n0Down.GT.0) THEN
          DO n=1, n0Down
             WRITE(6,'(4a14,3g14.6,3a14,g14.6,a14)') &
                '?', '?', '?', '?', xp(1:3,iDown(n)), '?', '?', '?', sDown(n), '?'
          END DO
  END IF
  IF (n0UP.GT.1) THEN
          DO n=1, n0Up
             WRITE(6,'(7a14,3g14.6,a14,g14.6)') &
                '?', '?', '?', '?', '?', '?', '?', xp(1:3,iUp(n)), '?', sUp(n)
          END DO
  END IF

  duOld=0.d0
  DO n=1, Min(nUp-n0Up, nDown-n0Down)
     s = 0.5d0*( sUp(n+n0Up) + sDown(n+n0Down) )
     du(1:3) = u(1:3,iUp(n+n0Up)) - u(1:3,iDown(n+n0Down))
     ! Apply periodic boundary conditions 
     DO i=1, 3
        IF (uPerio(i).LE.distance_zero) Cycle
        IF (n.EQ.1) THEN ! First point
                !du(i) = Modulo( du(i), uPerio(i) )
                du(i) = du(i) - aNInt(du(i)/uPerio(i))*uPerio(i)
        ELSE
                DO 
                   IF (Abs(du(i)-uPerio(i)-duOld(i)).LT.Abs(du(i)-duOld(i))) THEN
                           du(i) = du(i) - uPerio(i)
                   ELSE IF (Abs(du(i)+uPerio(i)-duOld(i)).LT.Abs(du(i)-duOld(i))) THEN
                           du(i) = du(i) + uPerio(i)
                   ELSE
                           exit
                   END IF
                END DO
        END IF
     END DO
     duOld(:) = du(:)
     WRITE(6,'(12g14.6)') s, du(1:3), xp(1:3,iDown(n+n0Down)), xp(1:3,iUp(n+n0Up)), &
        sDown(n+n0Down), sUp(n+n0Up)
  END DO

  IF (nDown-n0Down.GT.nUp-n0Up) THEN
          DO n=nUp-n0Up+n0Down+1, nDown
             WRITE(6,'(4a14,3g14.6,3a14,g14.6,a14)') &
                '?', '?', '?', '?', xp(1:3,iDown(n)), '?', '?', '?', sDown(n), '?'
          END DO
  ELSE IF (nUP-n0up.GT.nDown-n0Down) THEN
          DO n=nDown-n0Down+n0Up+1, nUp
             WRITE(6,'(7a14,3g14.6,a14,g14.6)') &
                '?', '?', '?', '?', '?', '?', '?', xp(1:3,iUp(n)), '?', sUp(n)
          END DO
  END IF


  ! Calculate average displacement on upper and lower planes
  uDown(:) = 0.d0
  DO n=1, nDown
     uDown(:) = uDown(:) + u(1:3,iDown(n))
  ENDDO
  uUp(:) = 0.d0
  DO n=1, nUp
     uUp(:) = uUp(:) + u(1:3,iUp(n))
  ENDDO

  WRITE(6,*)
  WRITE(6,'(2(a,i0))') "# nUp = ", nUp, '  -  nDown = ', nDown
  WRITE(6,'(a,3g14.6)') "# dU = ", uUp(:) - uDown(:)

  DEALLOCATE(sUp) ; DEALLOCATE(sDown)
  DEALLOCATE(iUp) ; DEALLOCATE(iDown)
  DEALLOCATE(xp) ; DEALLOCATE(u)

CONTAINS

  FUNCTION CrossProduct(ex,ey) RESULT(ez)
    
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: ex, ey
    REAL(kind(0.d0)), dimension(1:3) :: ez

    ez(1) = ex(2)*ey(3) - ex(3)*ey(2)
    ez(2) = ex(3)*ey(1) - ex(1)*ey(3)
    ez(3) = ex(1)*ey(2) - ex(2)*ey(1)

  END FUNCTION CrossProduct

END PROGRAM Burgers1D
