MODULE Rearrange

  USE LineCouple_elasticity_ani, ONLY : max_nlc  
  USE Disloc_elasticity_ani, ONLY : max_nd

  SAVE

  !===============================================================================
  INTEGER :: nds        ! Number of different dislocation systems
  INTEGER, parameter :: max_nds=max_nd ! Maximal value for nds 
  REAL(kind(0.d0)), dimension(1:3,1:max_nds) :: lDisloSystem    ! Line vector
  ! Number of dislocations in each system
  INTEGER, dimension(1:max_nds) :: nDisloSystem
  ! Indexes of all dislocation belonging to a given system
  INTEGER, dimension(1:max_nd,1:max_nds) :: iDisloSystem
  ! System index which each cluster belongs to
  INTEGER, dimension(1:max_nd) :: disloSystem   
  ! bClosedSystem(m)=.true. if total Burgers vector is 0
  ! bDipoleSystem(m)=.true. if dislocations can be grouped in dipoles
  LOGICAL, dimension(1:max_nds) :: bClosedSystem, bDipoleSystem
  ! iDisloDipole(n1) = n2 if dislocations n1 and n2 belong to the same dipole
  !                    0 if the dipole cannot be defined
  INTEGER, dimension(1:max_nd) :: iDisloDipole

  !===============================================================================
  INTEGER :: nlcs        ! Number of different systems for line-force couple
  INTEGER, parameter :: max_nlcs=max_nlc  ! Maximal value for nlcs
  REAL(kind(0.d0)), dimension(1:3,1:max_nlcs) :: lLineCoupleSystem    ! Line vector
  ! Number of line-force couples in each system
  INTEGER, dimension(1:max_nlcs) :: nLineCoupleSystem
  ! Indexes of all line-force couples belonging to a given system
  INTEGER, dimension(1:max_nlc,1:max_nlcs) :: iLineCoupleSystem
  ! System index which each cluster belongs to
  INTEGER, dimension(1:max_nlc) :: LineCoupleSystem   

  !===============================================================================
  LOGICAL, dimension(1:max_nd, 1:max_nlc) :: dislo_lineCouple
  !===============================================================================
  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE RearrangeDislo(out)

    USE babel_data, ONLY : xImages, yImages, zImages, at, &
        verbosity, distance_zero2, distance_zero, matid
    USE Disloc_elasticity_ani
    USE Math
    IMPLICIT NONE

    INTEGER, intent(in), optional :: out

    INTEGER :: n, m, n1, n2, n3, np, i
    REAL(kind(0.d0)) :: inv_norm, scalar, ex_norm2
    REAL(kind(0.d0)), dimension(1:3) :: ex, ey, ez
    REAL(kind(0.d0)), dimension(1:3,1:max_nds) :: bDisloSystem
    REAL(kind(0.d0)), dimension(1:3) :: b1, b2, dA
    LOGICAL :: create
    CHARACTER(len=20) :: for

    lDisloSystem(1:3,:)=0.d0
    nds=0
    nDisloSystem(:)=0
    iDisloSystem(:,:)=0
    disloSystem(:)=0
    iDisloDipole(:)=0

    IF (nd.LE.0) Return

    IF ( (Present(out)) .AND.(verbosity.GE.verbosity_max) ) THEN
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
            WRITE(out,'(a)') 'Rearrange dislo in systems sharing the same line vector'
            WRITE(out,*)
    END IF
    
    dislo_loop1: DO n=1, nd

       !------------------------------------
       ! Normalize dislo line vector 
       inv_norm = 1.d0/sqrt( Sum( lDislo(1:3,n)**2 ) )
       lDislo(1:3,n) = inv_norm*lDislo(1:3,n)

       !------------------------------------
       ! Check cutting direction
       IF ( Abs(Sum(cutDislo(1:3,n)*lDislo(1:3,n))).GT.Distance_Zero ) THEN
               WRITE(0,'(a,i0)') 'Dislocation ', n
               WRITE(0,'(a)') 'Cutting direction cutDislo(1:3) has to be &
                   &orthogonal to the dislocation line vector'
               WRITE(0,'(a,3g14.6)') '  cutDislo(1:3) = ', cutDislo(1:3,n)
               WRITE(0,'(a,3g14.6)') '  lDislo(1:3) = ', lDislo(1:3,n)
               STOP '< RearrangeDislo >'
       END IF
       ex_norm2 = SUM( cutDislo(1:3,n)**2 )
       IF (ex_norm2.LE.Distance_Zero2) THEN
               ! Try edge component of Burgers vector
               ex(1:3) = bDislo(1:3,n) -  Sum( bDislo(1:3,n)*lDislo(1:3,n) ) * lDislo(1:3,n)
               ex_norm2 = Sum( ex(1:3)**2 )
               IF (ex_norm2.LE.Distance_Zero2) THEN
                       ! Try vector perpendicular to ldislo(1:3) and (1,0,0)
                       ex(1)=0.d0 ; ex(2)=lDislo(3,n) ; ex(3)=-lDislo(2,n)
                       ex_norm2 = Sum( ex(1:3)**2 )
                       IF (ex_norm2.LE.Distance_Zero2) THEN
                               ! Try vector perpendicular to ldislo(1:3) and (0,1,0)
                               ex(1)=-lDislo(3,n) ; ex(2)=0.d0 ; ex(3)=lDislo(1,n)
                               ex_norm2 = Sum( ex(1:3)**2 )
                       ENDIF
               END IF
               cutDislo(1:3,n) = ex(1:3)/Sqrt(ex_norm2)
               IF (Present(out)) THEN
                       WRITE(out,*)
                       WRITE(out,'(a,i0)') 'Dislocation ', n
                       WRITE(out,'(a)') '  cutting direction has been initialized'
                       WRITE(out,'(a,3g14.6)') '  cutDislo(1:3) = ', cutDislo(1:3,n)
               END IF
       END IF

       !------------------------------------
       ! Organize dislo in different systems
       ! a system is defined by a common line direction
       create=.true.
       system_loop1: DO m=1, nds
          ! Scalar product between the line directions of the system and of the
          ! dislo
          scalar = Sum( lDisloSystem(1:3,m)*lDislo(1:3,n) )
          IF (Abs(scalar-1.d0).LE.distance_zero) THEN
                  ! Dislo n belongs to system m
                  disloSystem(n)=m
                  nDisloSystem(m) = nDisloSystem(m)+1
                  iDisloSystem(nDisloSystem(m),m) = n
                  create=.false.
                  CYCLE
          ELSE IF (Abs(scalar+1.d0).LE.distance_zero) THEN
                  ! Dislo n belongs to system m
                  ! We need to invert its orientation
                  disloSystem(n)=m
                  nDisloSystem(m) = nDisloSystem(m)+1
                  iDisloSystem(nDisloSystem(m),m) = n
                  lDislo(1:3,n) = -lDislo(1:3,n)
                  bDislo(1:3,n) = -bDislo(1:3,n)
                  create=.false.
                  CYCLE
          END IF
       END DO system_loop1
       IF (create) THEN
               ! We need to create a new dislo system
               nds=nds+1
               disloSystem(n)=nds
               lDisloSystem(1:3,nds) = lDislo(1:3,n)
               nDisloSystem(m) = 1
               iDisloSystem(1,m) = n
       END IF

    END DO dislo_loop1

    !------------------------------------
    ! Check if total Burgers vector of a system equals 0.
    ! If this is the case, try to group dislocations as dipoles
    system_loop2: DO m=1, nds ! Loop on all systems
       bDisloSystem(1:3,m) = 0.d0
       DO n=1, nDisloSystem(m) ! Loop on all dislocations belonging to the system
          ! Sum of Burgers vectors
          bDisloSystem(1:3,m) = bDisloSystem(1:3,m) + bDislo(1:3,iDisloSystem(n,m))
       END DO
       bClosedSystem(m) = ( Sum( bDisloSystem(1:3,m)**2 ) .LE.distance_Zero2)
       bDipoleSystem(m) = .false.
       IF (bClosedSystem(m)) THEN       ! Try to rearrange dislo by dipole pair
               IF (nDisloSystem(m).EQ.2) THEN
                       bDipoleSystem(m)=.true.
               ELSE
                       DO n=1, nDisloSystem(m), 2
                          bDipoleSystem(m)=.FALSE. ! Initialization
                          ! Index and Burgers vector of dislo 1
                          n1 = iDisloSystem(n,m)
                          b1(:) = bDislo(1:3,n1)
                          DO np=n+1, nDisloSystem(m)
                             n2 = iDisloSystem(np,m)
                             b2(:) = b1(:) + bDislo(1:3,n2)
                             IF ( Sum( b2(:)**2 ) .LE.distance_Zero2) THEN
                                     ! Dislo n1 and n2 form a dipole
                                     IF (np.NE.n+1) THEN
                                             n3 = iDisloSystem(n+1,m)
                                             iDisloSystem(n+1,m) = iDisloSystem(np,m)
                                             iDisloSystem(np,m) = n3
                                     END IF
                                     bDipoleSystem(m)=.TRUE.
                                     iDisloDipole(n1)=n2
                                     iDisloDipole(n2)=n1
                                     EXIT
                             END IF
                          END DO
                          IF (.NOT.bDipoleSystem(m)) EXIT
                       END DO
               END IF
       END IF

       ! Check compatibility of cutting directions with dipole definitions
       IF ( bDipoleSystem(m)) THEN
               DO n=1, nDisloSystem(m), 2
                  n1 = iDisloSystem(n,m)
                  n2 = iDisloSystem(n+1,m)
                  dA(:) = cDislo(:,n1) - cDislo(:,n2)
                  IF ( .NOT.Colinear(dA(:),cutDislo(:,n1)) ) THEN
                                  cutDislo(:,n1) = dA(:)
                                  IF ( (Present(out)).AND.(verbosity.GE.verbosity_max) ) THEN
                                          WRITE(out,'(3(a,i0))') '  Cutting direction of dislocation ', &
                                                n1, ' modified to correspond to dipole vector ', &
                                                n1, ' - ', n2
                                  END IF
                  END IF
                  IF ( .NOT.Colinear(dA(:),cutDislo(:,n2)) ) THEN
                                  cutDislo(:,n2) = dA(:)
                                  IF ( (Present(out)).AND.(verbosity.GE.verbosity_max) ) THEN
                                          WRITE(out,'(3(a,i0))') '  Cutting direction of dislocation ', &
                                                n2, ' modified to correspond to dipole vector ', &
                                                n1, ' - ', n2
                                  END IF
                  END IF
                  ! Periodic boundary conditions and the calculation of the
                  ! homogeneous strain do not work if the dipole cut
                  ! is not contained inside the simulation box
                  scalar = Sum( cutDislo(1:3,n1)*cutDislo(1:3,n2) )
                  IF (scalar.LT.0.d0) THEN
                          cutDislo(:,n2) = -cutDislo(:,n2)
                          IF ( (Present(out)).AND.(verbosity.GE.verbosity_max) ) THEN
                                  WRITE(out,'(a,i0,a)') '  Cutting direction of dislocation ', &
                                        n2, ' reversed'
                          END IF
                         !!$WRITE(0,*)
                         !!$WRITE(0,'(a)') 'WARNING !!!!'
                         !!$WRITE(0,'(2(a,i0),a)') 'Dislocations ', n1, ' and ', n2, &
                                !!$' have a cut surface not contained in the primitive unit cell'
                         !!$WRITE(0,'(a)') 'Periodic boundary conditions may not work'
                         !!$WRITE(0,'(a,i0,a,3g14.6)') '  cutDislo(1:3,',n1,') = ', cutDislo(1:3,n1)
                         !!$WRITE(0,'(a,i0,a,3g14.6)') '  cutDislo(1:3,',n2,') = ', cutDislo(1:3,n2)
                         !!$WRITE(0,*)
                         !!$WRITE(0,'(a)')'Use it at your own risks !!!!!!'
                         !!$WRITE(0,'(a)') '  (it may be necessary for grain boundaries)'
                         !!$DO i=5, 1, -1
                            !!$WRITE(0,*) i
                            !!$CALL Sleep(1)
                         !!$END DO
                  END IF
               END DO
       END IF
    END DO system_loop2

    !------------------------------------
    dislo_loop2: DO n=1, nd
 
       ! Make dislocation cut and line directions orthogonal
       scalar = Sum( cutDislo(:,n)*lDislo(:,n) )
       cutDislo(:,n) = cutDislo(:,n) - scalar*lDislo(:,n)

       ! Normalize dislocation cut direction
       inv_norm = 1.d0/Sqrt( SUM( cutDislo(1:3,n)**2 ) )
       cutDislo(1:3,n) = inv_norm*cutDislo(1:3,n)

       ! Orientate the crystal
       !  z: dislo vector line
       !  y: glide plane normal
       !  x: cutting direction corresponding to the discontinuity
       ez(1:3) = lDislo(1:3,n)
       ex(1:3) = cutDislo(1:3,n)
       ey(1:3) = CrossProduct( ez(1:3), ex(1:3) ) 

       ! Matrix to change the orientation of the crystal from cartesian coordinates
       ! (cubic axes) to the ref. frame of the dislocation
       rotDislo(1,1:3,n) = ex(1:3)
       rotDislo(2,1:3,n) = ey(1:3)
       rotDislo(3,1:3,n) = ez(1:3)

       ! Inverse matrix
       inv_rotDislo(1:3,1:3,n) = Transpose( rotDislo(1:3,1:3,n) )
       
       ! Check that rotDislo is a rotation matrix
       IF (MatNorm2( MatMul( inv_rotDislo(:,:,n), rotDislo(:,:,n) ) - matId ).GT.1.d-6) THEN
               WRITE(0,'(a,i0,a)') 'Rotation matrix for dislo ', n, ' is not a unitary matrix'
               STOP '< Rearrange >'
       END IF

    END DO dislo_loop2

    !CALL Init_iDisloDipole
    

    !------------------------------------
    IF ( (Present(out)).AND.(verbosity.GE.verbosity_max) ) THEN
            WRITE(out,*)
            WRITE(out,'(a,i0)') '  Number of different systems: ', nds
            IF (nds.GT.1) THEN
                    WRITE(out,'(a)') '    more than 1 system => no energy calculation'
            END IF
            DO m=1, nds
               WRITE(for,'(a,i0,a)') '(a,i0,a,', nDisloSystem(m), '(i0,1x))'
               WRITE(out,for) '    Dislocations belonging to system ', m, &
                        ': ', iDisloSystem(1:nDisloSystem(m),m)
               WRITE(out,'(a,3g14.6)') '      corresponding line vector: ', &
                        lDisloSystem(1:3,m)
               IF (bClosedSystem(m)) THEN
                       WRITE(out,'(a)') '      total Burgers vector is null'
                       IF (bDipoleSystem(m)) THEN
                              WRITE(out,'(a)') '      dislocations grouped in dipoles allowing energy calculation'
                       ELSE
                              WRITE(out,'(a)') '      dislocations not grouped in dipoles  &
                                        &(no energy calculation)'
                       END IF
               ELSE
                       WRITE(out,'(a,3(g14.6),a)') &
                                '      non-null total Burgers vector: b(1:3) = ', &
                                bDisloSystem(1:3,m), '  (no energy calculation)'
               END IF
            END DO
    END IF

    !------------------------------------
    ! Check that periodic boundary conditions are compatible with 
    ! line directions and Burgers circuits for all systems
    IF (xImages.OR.yImages.OR.zImages) THEN
            DO m=1, nds ! Loop on all systems
               IF ( ( xImages .AND. ( Colinear(at(:,1),lDisloSystem(:,m)) ) ).OR. &
                    ( yImages .AND. ( Colinear(at(:,2),lDisloSystem(:,m)) ) ).OR. &
                    ( zImages .AND. ( Colinear(at(:,3),lDisloSystem(:,m)) ) ) ) THEN
                       WRITE(0,'(a,i0,a)') 'Line direction of dislocations belonging to system ', &
                                m, ' is colinear with a periodicity vector'
                       WRITE(0,'(a)') 'One of the following variables should be put to .FALSE.: '
                       WRITE(0,'(a)') '   xImages, yImages or zImages'
                       STOP '< RearrangeDislo >'
               END IF
               IF (.NOT.bClosedSystem(m)) THEN
                       WRITE(0,*)
                       WRITE(0,'(a)') 'WARNING'
                       WRITE(0,'(a,i0,a,3g14.6)') 'Total Burgers vector corresponding to&
                                & system ', m, ': ', bDisloSystem(1:3,m)
                       WRITE(0,'(a)') 'This should be zero to be fully compatible with&
                                & periodic boundary conditions'
                       WRITE(0,*)
                       WRITE(0,'(a)')'Use it at your own risks !!!!!!'
                       WRITE(0,'(a)') '  (it may be necessary for grain boundaries)'
                       DO i=5, 1, -1
                          WRITE(0,*) i
                          CALL Sleep(1)
                       END DO
                       !!$STOP '< RearrangeDislo >'
               END IF
            END DO
    END IF

  END SUBROUTINE RearrangeDislo

  SUBROUTINE Init_iDisloDipole
    ! Initialize table iDisloDipole
    !   iDisloDipole(n1) = n2 if dislocations n1 and n2 belong to the same dipole
    !                      0 if the dipole cannot be defined
    !  ==> NOT USED: table is directly initialized in subroutine RearrangeDislo

    USE Disloc_elasticity_ani, ONLY : nd
    IMPLICIT NONE

    INTEGER :: n1, n2, m, i1, i2

    iDisloDipole(:) = 0

    DO n1=1, nd

            ! Dislo n1 belongs to system m
            m = disloSystem(n1)

            ! Check if dislocations are grouped in dipoles for this system
            IF (.NOT.bDipoleSystem(m)) Cycle

            ! Search dislo index i1 in dipole system m for dislo n1
            DO i1=1, nDisloSystem(m)
               IF (iDisloSystem(i1,m).EQ.n1) EXIT
            END DO

            ! Search dislo index i2 in dipole system m for dislo n2
            IF (Modulo(i1,2).EQ.0) THEN
                    i2 = i1 - 1
            ELSE
                    i2 = i1 + 1
            END IF

            ! Corresponding dislocation number
            iDisloDipole(n1) = iDisloSystem(i2,m)

    END DO

  END SUBROUTINE Init_iDisloDipole

  SUBROUTINE RearrangeLineCouple(out)
  
    USE babel_data, ONLY : xImages, yImages, zImages, at, &
        verbosity, distance_zero, matid
    USE LineCouple_elasticity_ani
    USE Math
    IMPLICIT NONE

    INTEGER, intent(in) :: out

    INTEGER :: n, m
    REAL(kind(0.d0)) :: inv_norm, scalar, ex_norm2, ey_norm2
    REAL(kind(0.d0)), dimension(1:3) :: ex, ey, ez
    LOGICAL :: create
    CHARACTER(len=20) :: for

    lLineCoupleSystem(1:3,:)=0.d0
    nlcs=0
    nLineCoupleSystem(:)=0
    iLineCoupleSystem(:,:)=0
    LineCoupleSystem(:)=0

    IF (nlc.LE.0) Return

    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
            WRITE(out,'(a)') 'Rearrange line-force couples in systems sharing the same line vector'
            WRITE(out,*)
    END IF
    
    LineCouple_loop1: DO n=1, nlc

       !------------------------------------
       ! Normalize line-force couple line vector 
       inv_norm = 1.d0/sqrt( Sum( lLineCouple(1:3,n)**2 ) )
       lLineCouple(1:3,n) = inv_norm*lLineCouple(1:3,n)

       !------------------------------------
       ! Check that force moments and line direction are orthogonal
       IF ( ( Abs( Sum( xLineCouple(1:3,n)*lLineCouple(1:3,n) ) ).GT.Distance_Zero ) .OR. &
        ( Abs( Sum( yLineCouple(1:3,n)*lLineCouple(1:3,n) ) ).GT.Distance_Zero ) ) THEN
               WRITE(0,'(a,i0)') 'Line-force couple ', n
               WRITE(0,'(a)') 'Force moment directions have to be orthogonal &
                   &to the line direction'
               STOP '< RearrangeLineCouple >'
       END IF

       ! Check that both force moments are orthogonal
       IF ( Abs( Sum( xLineCouple(1:3,n)*yLineCouple(1:3,n) )).GT. Distance_Zero ) THEN
               WRITE(0,'(a,i0)') 'Line-force couple ', n
               WRITE(0,'(a)') 'Force moment directions have to be orthogonal &
                   &to each other'
               STOP '< RearrangeLineCouple >'
       END IF

       !------------------------------------
       ! Organize line-force couple in different systems
       ! a system is defined by a common line direction
       create=.true.
       system_loop1: DO m=1, nlcs
          ! Scalar product between the line directions of the system
          ! and of the line-force couple
          scalar = Sum( lLineCoupleSystem(1:3,m)*lLineCouple(1:3,n) )
          IF (Abs(scalar-1.d0).LT.distance_zero) THEN
                  ! Couple n belongs to system m
                  LineCoupleSystem(n)=m
                  nLineCoupleSystem(m) = nLineCoupleSystem(m)+1
                  iLineCoupleSystem(nLineCoupleSystem(m),m) = n
                  create=.false.
                  CYCLE
          ELSE IF (Abs(scalar+1.d0).LT.distance_zero) THEN
                  ! Line-force couple n belongs to system m
                  ! We need to invert its orientation
                  LineCoupleSystem(n)=m
                  nLineCoupleSystem(m) = nLineCoupleSystem(m)+1
                  iLineCoupleSystem(nLineCoupleSystem(m),m) = n
                  lLineCouple(1:3,n) = -lLineCouple(1:3,n)
                  create=.false.
                  CYCLE
          END IF
       END DO system_loop1
       IF (create) THEN
               ! We need to create a new system
               nlcs=nlcs+1
               LineCoupleSystem(n)=nlcs
               lLineCoupleSystem(1:3,nlcs) = lLineCouple(1:3,n)
               nLineCoupleSystem(m) = 1
               iLineCoupleSystem(1,m) = n
       END IF

    END DO LineCouple_loop1

    !------------------------------------
    LineCouple_loop2: DO n=1, nlc
       !------------------------------------
       ! Orientate the crystal
       !  z: line-force couple direction
       ! Check x and y directions 
       ez(1:3) = lLineCouple(1:3,n)
       ex_norm2 = Sum( xLineCouple(1:3,n)**2 )
       IF (ex_norm2.NE.0.d0) THEN
               ex(1:3) = xLineCouple(1:3,n)/Sqrt(ex_norm2)
               ey(1:3) = CrossProduct( ez(1:3), ex(1:3) )
       ELSE
               ey_norm2 = Sum( yLineCouple(1:3,n)**2 )
               IF (ey_norm2.NE.0.d0) THEN
                       ey(1:3) = yLineCouple(1:3,n)/Sqrt(ey_norm2)
                       ex(1:3) = CrossProduct( ey(1:3), ez(1:3) )
               ELSE
                       WRITE(0,'(a,i0)') 'Line-force couple ', n
                       WRITE(0,'(a)') 'Directions for both moments are null'
                       STOP '< RearrangeLineCouple >'
               END IF
       END IF
       ! Matrix to change the orientation of the crystal from cartesian coordinates
       ! (cubic axex) to the ref. frame of the line force couple
       rotLineCouple(1,1:3,n) = ex(1:3)
       rotLineCouple(2,1:3,n) = ey(1:3)
       rotLineCouple(3,1:3,n) = ez(1:3)

       ! Inverse matrix
       inv_rotLineCouple(1:3,1:3,n) = Transpose( rotLineCouple(1:3,1:3,n) )

       ! Check that rotlineCouple is a rotation matrix
       IF (MatNorm2( MatMul( inv_rotlineCouple(:,:,n), rotlineCouple(:,:,n) ) - matId ).GT.1.d-6) THEN
               WRITE(0,'(a,i0,a)') 'Rotation matrix for line couple ', n, ' is not a unitary matrix'
               STOP '< Rotate >'
       END IF
    END DO LineCouple_loop2
    

    !------------------------------------
    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,'(a,i0)') '  Number of different systems: ', nlcs
            IF (nlcs.GT.1) THEN
                    WRITE(out,'(a)') '    more than 1 system => no energy calculation'
            END IF
            DO m=1, nlcs
               WRITE(for,'(a,i0,a)') '(a,i0,a,', nLineCoupleSystem(m), '(i0,1x))'
               WRITE(out,for) '    Line-force couples belonging to system ', m, &
                        ': ', iLineCoupleSystem(1:nLineCoupleSystem(m),m)
               WRITE(out,'(a,3g14.6)') '      corresponding line vector: ', &
                        lLineCoupleSystem(1:3,m)
            END DO
    END IF

    !------------------------------------
    ! Check that periodic boundary conditions are compatible with 
    ! line directions 
    IF (xImages.OR.yImages.OR.zImages) THEN
            DO m=1, nlcs ! Loop on all systems
               IF ( ( xImages .AND. ( Colinear(at(:,1),lLineCoupleSystem(:,m)) ) ).OR. &
                    ( yImages .AND. ( Colinear(at(:,2),lLineCoupleSystem(:,m)) ) ).OR. &
                    ( zImages .AND. ( Colinear(at(:,3),lLineCoupleSystem(:,m)) ) ) ) THEN
                       WRITE(0,'(a,i0,a)') 'Line direction of line-force couples belonging to system ', &
                                m, ' is colinear with a periodicity vector'
                       WRITE(0,'(a)') 'One of the following variables should be put to .FALSE.: '
                       WRITE(0,'(a)') '   xImages, yImages or zImages'
                       STOP '< RearrangeLineCouple >'
               END IF
            END DO
    END IF

  END SUBROUTINE RearrangeLineCouple

  SUBROUTINE Cross_dislo_lineCouple(out)
    ! Check if dislocations and line-force couples can be considered as the same
    ! line defect: they can be described with the same line direction and the
    ! same center
 
    USE Babel_data
    USE LineCouple_elasticity_ani
    USE Disloc_elasticity_ani
    USE Math
    IMPLICIT NONE

    INTEGER, intent(in) :: out
    INTEGER :: n1, n2

    dislo_lineCouple(:,:) = .FALSE.
    IF ( (nlc.LE.0).OR.(nd.LE.0) ) Return
    dislo_loop: DO n1=1, nd
       lineCouple_loop: DO n2=1, nlc
          dislo_lineCouple(n1,n2) = Colinear(lDislo(:,n1),lLineCouple(:,n2)) &
                .AND.  Colinear(lDislo(:,n1), cLineCouple(:,n2)-cDislo(:,n1))
       END DO lineCouple_loop
    END DO dislo_loop

    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
            WRITE(out,'(a)') 'Associate dislocations and line-force couples'
            WRITE (out,'(a)') '  (line-defects sharing the same line direction and the same origin)'
            WRITE(out,*)
            DO n1=1, nd
               DO n2=1, nlc
                  IF (dislo_lineCouple(n1,n2)) WRITE(out,'(2(a,i0))') &
                        "  Dislocation ", n1, " associated with line-force couple ", n2
               END DO
            END DO
            WRITE(out,*)
    END IF

  END SUBROUTINE Cross_dislo_lineCouple
  
END MODULE Rearrange
