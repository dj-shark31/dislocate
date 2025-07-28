MODULE loopModule

  SAVE

  ! Definition of dislocation loops
  LOGICAL :: l_loop       ! True if loops are defined
  INTEGER :: nLoop      ! Number of loops
  INTEGER, parameter :: max_nLoop=20       ! Maximal value for nLoops
  INTEGER, parameter :: max_iLoop=100      ! Maximal number of points in a loop
  ! bLoop(1:3,n): Burgers vector of loop n
  REAL(kind(0.d0)), dimension(1:3,1:max_nLoop) :: bLoop 
  ! iLoop(n): number of points in loop n
  INTEGER, dimension(1:max_nLoop) :: iLoop      
  ! xLoop(1:3,i,n): coordinate of point i belonging to loop n
  !   The first (i=0) and last point (i=iLoop(n)) should be the same to close
  !   the loop.
  REAL(kind(0.d0)), dimension(1:3,0:max_iLoop, 1:max_nLoop) :: xLoop    
  ! cLoop(1:3,n): center of loop n (needed to define cut surface)
  REAL(kind(0.d0)), dimension(1:3, 1:max_nLoop) :: cLoop    

CONTAINS

  SUBROUTINE InitLoops()

    IMPLICIT NONE
    l_loop=.FALSE.
    nLoop=0
    bLoop(:,:)=0.d0
    iLoop(:)=0
    xLoop(:,:,:)=0.d0
    cLoop(:,:)=0.d0

  END SUBROUTINE InitLoops

  SUBROUTINE ReadLoops(inp)

    IMPLICIT NONE
    INTEGER, intent(in) :: inp

    REAL(kind(0.d0)), parameter :: Zero2=1.d-12
    REAL(kind(0.d0)) :: xLoop_noise
    REAL(kind(0.d0)), dimension(1:3,1:max_iLoop) :: s
    INTEGER :: n, iln

    NAMELIST /loops/ nLoop, bLoop, iLoop, xLoop, cLoop, xLoop_noise

    ! Initialize loop definition
    xLoop_noise=0.d0
    CALL Random_Seed
    CALL InitLoops

    ! Read loop definition
    READ(inp,nml=loops)
    l_loop=.TRUE.

    ! Check number of loops
    IF (nLoop.GT.max_nLoop) THEN
            WRITE(0,'(a,i0)') 'nLoop = ', nLoop
            WRITE(0,'(a,i0)') 'maximum value allowed (max_nLoop) = ', max_nLoop
            STOP '< ReadLoop >'
    ELSE IF (nLoop.LE.0) THEN
            CALL InitLoops()
            RETURN
    END IF

    DO n=1, nLoop

       ! Number of points in the loop
       iln = iLoop(n)

       ! Check number of points for each loop
       IF (iln.GT.max_iLoop) THEN
               WRITE(0,'(3(a,i0))') 'Number of points for loop ', n, ': iLoop(', n,') = ', iln
               WRITE(0,'(a,i0)') 'maximum value allowed (max_iLoop) = ', max_iLoop
               STOP '< ReadLoop >'
       ELSE IF (iln.LE.0) THEN
               WRITE(0,'(3(a,i0))') 'Number of points for loop ', n, ': iLoop(', n,') = ', iln
               STOP '< ReadLoop >'
       END IF

       ! Check Burgers vector for each loop
       IF ( Sum(bLoop(:,n)**2).LE.Zero2 ) THEN
               WRITE(0,'(2(a,i0),a,3g20.12)') 'Burgers vector for loop ', n, ': bLoop(:,', n, ') = ', bLoop(:,n)
               WRITE(0,'(a,g20.12)') 'Norm smaller than Zero = ', Sqrt(Zero2)
               STOP '< ReadLoop >'
       END IF

       ! Add some noise to positions of loop points
       IF (xLoop_noise.GT.0.d0) THEN
               CALL Random_Number(s(:,1:iln))
               xLoop(:,1:iln,n) = xLoop(:,1:iln,n) + xLoop_noise*(s(:,1:iln)-0.5d0)
       END IF
    
       ! Define loop center if not already defined
       IF ( Sum(cLoop(:,n)**2).LE.Zero2 ) THEN
               ! It has not been defined a priori => we take the gravity center
               cLoop(:,n) = Sum( xLoop(:,1:iln,n), 2 )/dble(iln)
       ELSE IF (xLoop_noise.GT.0.d0) THEN
               CALL Random_Number(s(:,1))
               cLoop(:,n) = cLoop(:,n) + xLoop_noise*(s(:,1)-0.5d0)
       END IF

       ! Close the loop
       xLoop(:,0,n) = xLoop(:,iln,n)

    END DO
  
  END SUBROUTINE ReadLoops

  SUBROUTINE PrintLoops(out)
    ! Print loop definition on output unit out

    IMPLICIT NONE
    INTEGER, intent(in) :: out ! Output unit
    INTEGER :: n, i

    IF (.NOT.l_loop) THEN
            WRITE(out,'(a)') ' No loop defined'
            WRITE(out,*)
            RETURN
    END IF

    WRITE(out,'(a,i0)') ' Number of loops defined: nLoop = ', nLoop
    DO n=1, nLoop
       WRITE(out,'(a,i0)') ' Loop ' , n
       WRITE(out,'(a,3g20.12)') '   Burgers vector: ', bLoop(:,n)
       WRITE(out,'(a,3g20.12)') '   loop center: ', cLoop(:,n)
       WRITE(out,'(a,i0)')      '   number of points: ', iLoop(n)
       WRITE(out,'(a)')         '   coordinates of each point:'
       DO i=1, iLoop(n)
          WRITE(out,'(a,3g20.12)') '     ', xLoop(:,i,n)
       END DO
    END DO
    WRITE(out,*)

  END SUBROUTINE PrintLoops

  SUBROUTINE ScaleLoops(a)
    ! Multiply distances by a for loop definition

    IMPLICIT NONE
    REAL(kind(0.d0)), intent(in) :: a

    bLoop(:,1:nLoop) = a*bLoop(:,1:nLoop)
    xLoop(:,:,1:nLoop) = a*xLoop(:,:,1:nLoop)
    cLoop(:,1:nLoop) = a*cLoop(:,1:nLoop)

  END SUBROUTINE ScaleLoops

  SUBROUTINE RotateLoops(rot, out, verbose)
    ! Rotate loops according to rotation matrix rot(1:3,1:3)
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n

    bLoop(:,1:nLoop) = MatMul( rot(:,:), bLoop(:,1:nLoop) )
    DO n=1, nLoop
            xLoop(:,1:iLoop(n),n) = MatMul( rot(:,:), xLoop(:,1:iLoop(n),n) )
    END DO
    cLoop(:,1:nLoop) = MatMul( rot(:,:), cLoop(:,1:nLoop) )

    IF (verbose) THEN
            WRITE(out,'(a)') '  New definition for loops after rotation'
            CALL PrintLoops(out)
    END IF

  END SUBROUTINE RotateLoops

  SUBROUTINE TranslateLoops(u, out, verbose)
    ! Add displacement u(1:3) to loop defininition
    ! Print new definition on output unit out if verbose=.true.

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: u
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    INTEGER :: n, i

    DO n=1, nLoop
       DO i=1, iLoop(n)
          xLoop(:,i,n) = xLoop(:,i,n) + u(:)
       END DO
       cLoop(:,n) = cLoop(:,n) + u(:)
    END DO

    IF (verbose) THEN
            WRITE(out,'(a)') '  new definition for loops after translation'
            CALL PrintLoops(out)
    END IF

  END SUBROUTINE TranslateLoops

  FUNCTION Loops_displacement(r, n) RESULT(u)
    ! Calculate displacement created by loop n in R(1:3)

    USE Elasticity_loops
    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: r
    INTEGER, intent(in) :: n
    REAL(kind(0.d0)), dimension(1:3) :: u
    

    INTEGER :: i
    REAL(kind(0.d0)) :: omega
    REAL(kind(0.d0)), dimension(1:3) :: xA, xB, xC, b

    ! Poisson coefficient
    REAL(kind(0.d0)), parameter :: nu=1.d0/3.d0

    u(:) = 0.d0

    IF (n.GT.nLoop) Return

    ! Loop center
    xC(:) = cLoop(:,n) - r(:)

    ! Burgers vector
    b(:) = bLoop(:,n)

    omega=0.d0
    DO i=1, iLoop(n)
       xA(:) = xLoop(:,i-1,n) - r(:)
       xB(:) = xLoop(:,i,n) - r(:)
       ! Part due to solid angle
       omega = omega + SolidAngle(xA, xB, xC)
       ! Part due to elasticity
       u(:) = u(:) + DisloSeg_displacement_iso(xA, xB, b(:), nu)
    END DO
    u(:) = u(:) + bLoop(:,n)*omega

  END FUNCTION Loops_displacement

END MODULE loopModule
