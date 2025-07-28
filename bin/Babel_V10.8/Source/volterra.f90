MODULE Volterra_Module

        PRIVATE :: RemoveDuplicatedAtoms

CONTAINS

  SUBROUTINE RemoveAddDisloCut(xp0, xp, iTyp, im, keep, uAdd, at, out, verbose)
     ! Remove or add atoms in the cut created by dislocation of vacancy / interstitial type

     USE Math
     USE Babel_Data, ONLY : xImages, yImages, zImages, debug, distance_zero2
     USE Rearrange
     USE Disloc_elasticity_ani
     IMPLICIT NONE

     ! Atoms position: xp(1:3,i) for atom i (xp0: initial, xp: strained)
     REAL(kind(0.d0)), dimension(:,:), intent(inout) :: xp0, xp
     ! Atom type
     INTEGER, dimension(:), intent(inout) :: iTyp
     ! Number of atoms
     INTEGER, intent(inout) :: im
     ! Boolean saying if atoms should be kept or not
     LOGICAL, dimension(:), intent(out) :: keep
     ! uAdd(1:3,i : plastic displacement to add to atom i because of Volterra cut
     REAL(kind(0.d0)), dimension(:,:), intent(out) :: uAdd
     ! Periodicity vectors
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
     ! Boolean to control verbose mode
     LOGICAL, intent(in) :: verbose


     INTEGER, intent(in) :: out

     INTEGER :: i, n, nRemove, nAdd, imNew, imInitial, imm
     REAL(kind(0.d0)) ::  area, h, volume, Va
     REAL(kind(0.d0)), dimension(1:3) :: be, origin, s, dipole, xpAtom, xp0Atom
     REAL(kind(0.d0)), dimension(1:3,1:3) :: Hmat, Gmat
     LOGICAL :: vType

     ! Periodic boundary conditions
     INTEGER :: nPerio, np
     REAL(kind(0.d0)), dimension(1:3,1:27) :: uPerio


     ! Maximal number of atoms allowed by array size
     imm=Size(xp0,2)
     IF ((imm.LT.im).OR.(Size(iTyp).LT.im).OR.(Size(keep).LT.im).OR.(Size(uAdd,2).LT.im)) THEN
             WRITE(0,'(a)') 'Problem with size of arrays on input'
             WRITE(0,'(a,i0)') '  imm = ', imm
             WRITE(0,'(a,i0)') '  im  = ', im 
             STOP '< RemoveAddDisloCut >'
     END IF
     IF ((Size(xp,2).NE.imm).OR.(Size(iTyp).NE.imm).OR.(Size(keep).NE.imm).OR.(Size(uAdd,2).NE.imm)) THEN
             WRITE(0,'(a)') 'Problem with size of arrays on input'
             STOP '< RemoveAddDisloCut >'
     END IF

     IF ((.NOT.l_dislo).OR.(nd.LT.1)) Return

     ! Initialization
     keep(:)=.TRUE.
     uAdd(:,:)=0.d0
    
     ! Store initial number of atom
     imNew=im
     imInitial=im

     ! Define vectors for periodic boundary conditions
     uPerio(:,:)=0.d0
     nPerio=1
     IF ( xImages ) THEN
             DO np=1, nPerio
                uPerio(:,nPerio+np)   = uPerio(:,np) - at(:,1)
                uPerio(:,2*nPerio+np) = uPerio(:,np) + at(:,1)
             END DO
             nPerio = 3*nPerio
     END IF
     IF ( yImages ) THEN
             DO np=1, nPerio
                uPerio(:,nPerio+np)   = uPerio(:,np) - at(:,2)
                uPerio(:,2*nPerio+np) = uPerio(:,np) + at(:,2)
             END DO
             nPerio = 3*nPerio
     END IF
     IF ( zImages ) THEN
             DO np=1, nPerio
                uPerio(:,nPerio+np)   = uPerio(:,np) - at(:,3)
                uPerio(:,2*nPerio+np) = uPerio(:,np) + at(:,3)
             END DO
             nPerio = 3*nPerio
     END IF
     IF ( DEBUG ) THEN
             WRITE(out,*)
             WRITE(out,'(a,i0)') 'Remove Volterra cut - number of periodicity vectors: ', nPerio
             DO np=1, nPerio
                WRITE(out,'(2x,3g20.6)') uPerio(1:3,np)
             END DO
             WRITE(out,*)
     END IF
     
     ! ----- Isolated Dislocations -------------------------
     IF (nds.LE.0.OR..NOT.ALL(bDipoleSystem(1:nds))) THEN
             ! Dislocations are not grouped in dipoles
             !   => we need to check each dislocation cut
             dislo_loop1: DO n=1, nd

               nRemove=0
               nAdd=0
             
               ! Check if Burgers vector is colinear with line direction
               !        screw dislo => no atom needs to be removed
               IF ( Colinear( bDislo(:,n), lDislo(:,n) ) ) Cycle

               ! Edge component of the Burgers vector
               be(:) = bDislo(:,n) - Sum( bDislo(:,n)*lDislo(:,n) )*lDislo(:,n)

               ! Check if edge component of the Burgers vector is colinear with cut direction
               !        cut = glide plane => no atom needs to be removed
               IF ( Colinear( be(:), cutDislo(:,n) ) ) Cycle

               ! Check if the cut is of I or V type
               IF ( ScalarTripleProduct( cutDislo(:,n), be(:), lDislo(:,n)).LT.0.d0 ) THEN
                       vType=.true.
                       IF (verbose) WRITE(out,'(a)') 'Volterra cut of V type'
               ELSE
                       vType=.FALSE.
                       IF (verbose) WRITE(out,'(a)') 'Volterra cut of I type'
               END IF

               ! Matrix to change reference frame: ex = cutDislo(:,n)
               !                                   ey = be(:)
               !                                   ez = lDislo(:,n)
               Hmat(:,1) = cutDislo(:,n)
               Hmat(:,2) = be(:)
               Hmat(:,3) = lDislo(:,n)
               CALL MatInv(Hmat,Gmat)
               Origin(:) = cDislo(:,n) 

               DO i=1, imNew ! Loop on atoms
                  xp0Atom(:) = xp0(:,i) 
                  xpAtom(:)  = xp(:,i)  
                  ! Atoms reduced coordinates in new reference frame
                  s(:) = MatMul( Gmat(:,:), xpAtom(:)-Origin(:) )
                  IF ( (s(1).LE.0.d0).AND.(s(2).GT.-0.5d0).AND.(s(2).LT.0.5d0) ) THEN
                          IF (vType) THEN
                                  ! Vacancy type
                                  IF (keep(i)) THEN
                                          keep(i)=.FALSE.
                                          nRemove = nRemove+1
                                  END IF
                          ELSE
                                  ! Interstitial type
                                  nAdd = nAdd+1
                                  im=im+1
                                  IF (im.GT.imm) THEN
                                          WRITE(0,'(a)') 'Size of atom arrays is too small'
                                          WRITE(0,'(a)') 'You need to increase value of imm in the input file'
                                          STOP '< RemoveAddDisloCut >'
                                  END IF
                                  ! Create new atom
                                  xp0(1:3,im) = xp0Atom(1:3)
                                  xp(1:3,im)  = xpAtom(1:3)
                                  iTyp(im) = iTyp(i)
                                  keep(im) = keep(i)
                                  IF (s(2).GT.0.d0) THEN 
                                          uAdd(:,im) = uAdd(:,im) - 1.d0*bDislo(:,n)
                                  ELSE                            
                                          uAdd(:,im) = uAdd(:,im) + 1.d0*bDislo(:,n)
                                  END IF
                                  IF (DEBUG) THEN
                                          WRITE(out,'(2(a,i0))') "Add atom ", i, " as new atom ", im
                                          WRITE(out,'(a,3g20.12)') "  position: xp0(1:3) = ", xp0(1:3,im)
                                          WRITE(out,'(a,3g20.12)') "  added displacement: uAdd(1:3) = ", uAdd(1:3,im)
                                  END IF
                          END IF
                  END IF

               END DO

             END DO dislo_loop1

     ! ----- Dislocation Dipoles -------------------------
     ELSE

             nRemove=0
             nAdd=0

             ! Dislocations are grouped in dipoles
             !   => we need to check each dipole cut
             dislo_loop2: DO n=1, nd, 2

               ! Check if Burgers vector is colinear with line direction
               !        screw dislo => no atom needs to be removed
               IF ( Colinear( bDislo(:,n), lDislo(:,n) ) ) Cycle

               ! Edge component of the Burgers vector
               be(:) = bDislo(:,n) - Sum( bDislo(:,n)*lDislo(:,n) )*lDislo(:,n)

               ! Vector joining centers of both dislocations composing the dipole
               dipole(:) = cDislo(:,n+1) - cDislo(:,n)

               ! Check if edge component of the Burgers vector is colinear with dipole direction
               !        cut = glide plane => no atom needs to be removed
               IF ( Colinear( be(:), dipole(:) ) ) Cycle

               ! Calculate volume that needs to be removed
               area =  ScalarTripleProduct( dipole(:), be(:), lDislo(:,n) )
               IF ( Abs(area).LE.distance_zero2) CYCLE
               IF ( Colinear( lDislo(:,n), at(:,1) ) ) THEN
                       h = Abs( Sum( lDislo(:,n)*at(:,1) ) )
               ELSE IF ( Colinear( lDislo(:,n), at(:,2) ) ) THEN
                       h = Abs( Sum( lDislo(:,n)*at(:,2) ) )
               ELSE IF ( Colinear( lDislo(:,n), at(:,3) ) ) THEN
                       h = Abs( Sum( lDislo(:,n)*at(:,3) ) )
               END IF
               volume = h*area

               ! Check if the cut is of I or V type
               IF ( area.GT.0.d0 ) THEN
                       vType=.true.
                       IF (verbose) WRITE(out,'(a)') 'Volterra cut of V type'
               ELSE
                       vType=.FALSE.
                       IF (verbose) WRITE(out,'(a)') 'Volterra cut of I type'
               END IF

               ! Matrix to change reference frame: ex = cutDislo(:,n)
               !                                   ey = be(:)
               !                                   ez = lDislo(:,n)
               Hmat(:,1) = dipole(:)
               Hmat(:,2) = be(:)
               Hmat(:,3) = lDislo(:,n)
               CALL MatInv(Hmat,Gmat)
               Origin(:) = cDislo(:,n) 

               atom_loop2: DO i=1, imNew
                  xp0Atom(:) = xp0(:,i) 
                  xpAtom(:)  = xp(:,i)  
                  perio_loop2: DO np=1, nPerio
                     ! Atoms reduced coordinates in new reference frame
                     s(:) = MatMul( Gmat(:,:), xpAtom(:)-Origin(:) + uPerio(:,np) )
                     IF ( (s(1).EQ.0.d0).OR.(s(1).EQ.1.d0).OR.(s(2).EQ.-0.5d0).OR.(s(2).EQ.0.5d0) ) THEN
                             WRITE(out,'(a)') 'One of the atom is located exactly at &
                                   &the frontier of the volume to remove in Volterra procedure'
                             WRITE(out,'(a)') 'You need to add some noise on dislocation positions'
                             STOP '< RemoveAddDisloCut >' 
                     ELSE IF ( (s(1).GT.0.d0).AND.(s(1).LT.1.d0).AND.(s(2).GT.-0.5d0).AND.(s(2).LT.0.5d0) ) THEN
                             IF (vType) THEN
                                     ! Vacancy type
                                     IF (keep(i)) THEN
                                             keep(i)=.FALSE.
                                             nRemove = nRemove+1
                                     END IF
                             ELSE
                                     ! Interstitial type
                                     nAdd = nAdd+1
                                     im=im+1
                                     IF (im.GT.imm) THEN
                                             WRITE(0,'(a)') 'Size of atom arrays is too small'
                                             WRITE(0,'(a)') 'You need to increase value of imm in the input file'
                                             STOP '< RemoveAddDisloCut >' 
                                     END IF
                                     ! Create new atom
                                     xp0(1:3,im) = xp0Atom(1:3)
                                     xp(1:3,im)  = xpAtom(1:3)
                                     iTyp(im) = iTyp(i)
                                     keep(im) = keep(i)
                                     IF (s(2).GT.0.d0) THEN 
                                             uAdd(:,im) = uAdd(:,i) - 1.d0*bDislo(:,n)
                                     ELSE                            
                                             uAdd(:,im) = uAdd(:,i) + 1.d0*bDislo(:,n)
                                     END IF
                                     IF (DEBUG) THEN
                                             WRITE(out,'(2(a,i0))') "Add atom ", i, " as new atom ", im
                                             WRITE(out,'(a,3g20.12)') "  position: xp0(1:3) = ", xp0(1:3,im)
                                             WRITE(out,'(a,3g20.12)') "  added displacement: uAdd(1:3) = ", uAdd(1:3,im)
                                     END IF
                             END IF
                     END IF

                  END DO perio_loop2

               END DO atom_loop2

               ! New number of atoms for next loop
               imNew = im

               IF (verbose) THEN
                       WRITE(out,*)
                       WRITE(out,'(2(a,i0))') 'Volterra cut for dipole joining dislocation n=', n,&
                                ' and n+1=', n+1
                       WRITE(out,'(2x,2(i0,a))') nRemove, ' atoms removed  -  ', nAdd, ' atoms added'
                       ! Corresponding number of atoms that will be theoretically removed
                       Va = Abs( MatDet( at(:,:) ) )/dble(im)
                       WRITE(out,'(a,2(f0.2,a))') '  Volume of the Volterra cut, V = area * h, with area = ', area, &
                                ', and h = ', h 
                       WRITE(out,'(a,3(f0.2,a))') '                              V = ', area*h, &
                                ', corresponds to ', area*h/VA , ' atomic volumes ( Va = ', Va,' )'
                       WRITE(out,*)
               END IF

             END DO dislo_loop2
     END IF

     ! Remove atoms which have been inserted multiple times by dislocation branch cuts
     CALL RemoveDuplicatedAtoms(xp, im, keep, uAdd, out, verbose) 


  END SUBROUTINE RemoveAddDisloCut

  SUBROUTINE RemoveAddDDipoleCut(xp0, xp, iTyp, im, keep, uAdd, at, out, verbose)
     ! Remove or add atoms in the cut created by dislocation dipole of vacancy / interstitial type

     USE Math
     USE Babel_Data, ONLY : xImages, yImages, zImages, debug, distance_zero2
     USE DDipoleModule
     IMPLICIT NONE

     ! Atoms position: xp(1:3,i) for atom i (xp0: initial, xp: strained)
     REAL(kind(0.d0)), dimension(:,:), intent(inout) :: xp0, xp
     ! Atom type
     INTEGER, dimension(:), intent(inout) :: iTyp
     ! Number of atoms
     INTEGER, intent(inout) :: im
     ! Boolean saying if atoms should be kept or not
     LOGICAL, dimension(:), intent(out) :: keep
     ! uAdd(1:3,i : plastic displacement to add to atom i because of Volterra cut
     REAL(kind(0.d0)), dimension(:,:), intent(out) :: uAdd
     ! Periodicity vectors
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
     ! Boolean to control verbose mode
     LOGICAL, intent(in) :: verbose


     INTEGER, intent(in) :: out

     INTEGER :: i, n, nRemove, nAdd, imInitial, imNew, imm
     REAL(kind(0.d0)) ::  area, h, volume, Va
     REAL(kind(0.d0)), dimension(1:3) :: be, origin, s, dipole, xpAtom, xp0Atom
     REAL(kind(0.d0)), dimension(1:3,1:3) :: Hmat, Gmat
     LOGICAL :: vType

     ! Periodic boundary conditions
     INTEGER :: nPerio, np
     REAL(kind(0.d0)), dimension(1:3,1:27) :: uPerio

     ! Maximal number of atoms allowed by array size
     imm=Size(xp0,2)
     IF ((imm.LT.im).OR.(Size(iTyp).LT.im).OR.(Size(keep).LT.im).OR.(Size(uAdd,2).LT.im)) THEN
             WRITE(0,'(a)') 'Problem with size of arrays on input'
             WRITE(0,'(a,i0)') '  imm = ', imm
             WRITE(0,'(a,i0)') '  im  = ', im 
             STOP '< RemoveAddDDipoleCut >'
     END IF
     IF ((Size(xp,2).NE.imm).OR.(Size(iTyp).NE.imm).OR.(Size(keep).NE.imm).OR.(Size(uAdd,2).NE.imm)) THEN
             WRITE(0,'(a)') 'Problem with size of arrays on input'
             STOP '< RemoveAddDDipoleCut >'
     END IF

     IF ((.NOT.l_DDipole).OR.(nDDipole.LT.1)) Return

     ! Initialization
     keep(:)=.TRUE.
     uAdd(:,:)=0.d0
    
     ! Store initial number of atom
     imInitial=im
     imNew=im

     ! Define vectors for periodic boundary conditions
     uPerio(:,:)=0.d0
     nPerio=1
     IF ( xImages ) THEN
             DO np=1, nPerio
                uPerio(:,nPerio+np)   = uPerio(:,np) - at(:,1)
                uPerio(:,2*nPerio+np) = uPerio(:,np) + at(:,1)
             END DO
             nPerio = 3*nPerio
     END IF
     IF ( yImages ) THEN
             DO np=1, nPerio
                uPerio(:,nPerio+np)   = uPerio(:,np) - at(:,2)
                uPerio(:,2*nPerio+np) = uPerio(:,np) + at(:,2)
             END DO
             nPerio = 3*nPerio
     END IF
     IF ( zImages ) THEN
             DO np=1, nPerio
                uPerio(:,nPerio+np)   = uPerio(:,np) - at(:,3)
                uPerio(:,2*nPerio+np) = uPerio(:,np) + at(:,3)
             END DO
             nPerio = 3*nPerio
     END IF
     IF ( DEBUG ) THEN
             WRITE(out,*)
             WRITE(out,'(a,i0)') 'Remove Volterra cut - number of periodicity vectors: ', nPerio
             DO np=1, nPerio
                WRITE(out,'(2x,3g20.6)') uPerio(1:3,np)
             END DO
             WRITE(out,*)
     END IF
     
     ! ----- Dislocation Dipoles -------------------------

     ! Dislocations are grouped in dipoles
     !   => we need to check each dipole cut
     DO n=1, nDDipole

       nRemove=0
       nAdd=0

       ! Check if Burgers vector is colinear with line direction
       !        screw dislo => no atom needs to be removed
       IF ( Colinear( bDDipole(:,n), lDDipole(:,n) ) ) Cycle

       ! Edge component of the Burgers vector
       be(:) = bDDipole(:,n) - Sum( bDDipole(:,n)*lDDipole(:,n) )*lDDipole(:,n)

       ! Vector joining centers of both dislocations composing the dipole
       dipole(:) = c2DDipole(:,n) - c1DDipole(:,n)

       ! Check if edge component of the Burgers vector is colinear with dipole direction
       !        cut = glide plane => no atom needs to be removed
       IF ( Colinear( be(:), dipole(:) ) ) Cycle

       ! Calculate volume that needs to be removed
       area = ScalarTripleProduct( dipole(:), be(:), lDDipole(:,n) )
       IF ( Abs(area).LE.distance_zero2) CYCLE
       IF ( Colinear( lDDipole(:,n), at(:,1) ) ) THEN
               h = Abs( Sum( lDDipole(:,n)*at(:,1) ) )
       ELSE IF ( Colinear( lDDipole(:,n), at(:,2) ) ) THEN
               h = Abs( Sum( lDDipole(:,n)*at(:,2) ) )
       ELSE IF ( Colinear( lDDipole(:,n), at(:,3) ) ) THEN
               h = Abs( Sum( lDDipole(:,n)*at(:,3) ) )
       END IF
       volume = h*area

       ! Check if the cut is of I or V type
       IF ( area.GT.0.d0 ) THEN
               vType=.true.
               IF (verbose) WRITE(out,'(a)') 'Volterra cut of V type'
       ELSE
               vType=.FALSE.
               IF (verbose) WRITE(out,'(a)') 'Volterra cut of I type'
       END IF

       ! Matrix to change reference frame: ex = cutDislo(:,n)
       !                                   ey = be(:)
       !                                   ez = lDislo(:,n)
       Hmat(:,1) = dipole(:)
       Hmat(:,2) = be(:)
       Hmat(:,3) = lDDipole(:,n)
       CALL MatInv(Hmat,Gmat)
       Origin(:) = c1DDipole(:,n) 

       atom_loop: DO i=1, imNew ! Loop on atoms
          xp0Atom(:) = xp0(:,i) 
          xpAtom(:)  = xp(:,i)  
          perio_loop: DO np=1, nPerio
             ! Atoms reduced coordinates in new reference frame
             s(:) = MatMul( Gmat(:,:), xpAtom(:)-Origin(:) + uPerio(:,np) )
             IF ( (s(1).EQ.0.d0).OR.(s(1).EQ.1.d0).OR.(s(2).EQ.-0.5d0).OR.(s(2).EQ.0.5d0) ) THEN
                     WRITE(out,'(a)') 'One of the atom is located exactly at &
                           &the frontier of the volume to remove in Volterra procedure'
                     WRITE(out,'(a)') 'You need to add some noise on dislocation positions'
                     STOP '< RemoveAddDDipoleCut >' 
             ELSE IF ( (s(1).GT.0.d0).AND.(s(1).LT.1.d0).AND.(s(2).GT.-0.5d0).AND.(s(2).LT.0.5d0) ) THEN
                     IF (vType) THEN
                             ! Vacancy type
                             IF (keep(i)) THEN
                                     keep(i)=.FALSE.
                                     nRemove = nRemove+1
                             END IF
                     ELSE
                             ! Interstitial type
                             nAdd = nAdd+1
                             im=im+1
                             IF (im.GT.imm) THEN
                                     WRITE(0,'(a)') 'Size of atom arrays is too small'
                                     WRITE(0,'(a)') 'You need to increase value of imm in the input file'
                                     STOP '< RemoveAddDDipoleCut >' 
                             END IF
                             ! Create new atom
                             xp0(1:3,im) = xp0(1:3,i)
                             xp(1:3,im)  = xp(1:3,i)
                             iTyp(im) = iTyp(i)
                             keep(im) = keep(i)
                             IF (s(2).GT.0.d0) THEN 
                                     uAdd(:,im) = uAdd(:,i) - 1.d0*bDDipole(:,n)
                             ELSE                            
                                     uAdd(:,im) = uAdd(:,i) + 1.d0*bDDipole(:,n)
                             END IF
                             IF (DEBUG) THEN
                                     WRITE(out,'(2(a,i0))') "Add atom ", i, " as new atom ", im
                                     WRITE(out,'(a,3g20.12)') "  position: xp0(1:3) = ", xp0(1:3,im)
                                     WRITE(out,'(a,3g20.12)') "  added displacement: uAdd(1:3) = ", uAdd(1:3,im)
                                     WRITE(out,'(a,2(i0,1x),6g20.12)') '#VOLTERRA ', i, im, xp0(1:3,im), uAdd(1:3,im)
                             END IF
                     END IF
             END IF

          END DO perio_loop

       END DO atom_loop

       ! New number of atoms
       imNew = im

       IF (verbose) THEN
               WRITE(out,*)
               WRITE(out,'(a,i0)') 'Volterra cut for dislocation dipole n=', n
               WRITE(out,'(2x,2(i0,a))') nRemove, ' atoms removed  -  ', nAdd, ' atoms added'
               ! Corresponding number of atoms that will be theoretically removed
               Va = Abs( MatDet( at(:,:) ) )/dble(im)
               WRITE(out,'(a,2(f0.2,a))') '  Volume of the Volterra cut, V = area * h, with area = ', area, &
                        ', and h = ', h 
               WRITE(out,'(a,3(f0.2,a))') '                              V = ', area*h, &
                        ', corresponds to ', area*h/VA , ' atomic volumes ( Va = ', Va,' )'
               WRITE(out,*)
       END IF

     END DO

     ! Remove atoms which have been inserted multiple times by dislocation branch cuts
     CALL RemoveDuplicatedAtoms(xp, im, keep, uAdd, out, verbose) 

  END SUBROUTINE RemoveAddDDipoleCut

  SUBROUTINE RemoveAddLoopCut(xp0, xp, iTyp, im, keep, uAdd, at, out, verbose)
    ! Remove or add atoms in the loop cut for vacancy / interstitial prismatic loops

    USE Babel_data, ONLY : nxImages, nyImages, nzImages
    USE Math
    USE LoopModule
    IMPLICIT NONE
    
    LOGICAL, parameter :: Debug=.false.
    LOGICAL, parameter :: DebugPerio=.false.

    REAL(kind(0.d0)), parameter :: Zero_angle=1.d-6     ! radians
    REAL(kind(0.d0)), parameter :: Zero_angle2=Zero_angle*Zero_angle

    ! Atoms position: xp(1:3,i) for atom i (xp0: initial, xp: strained)
    REAL(kind(0.d0)), dimension(:,:), intent(inout) :: xp0, xp
    ! Atom type
    INTEGER, dimension(:), intent(inout) :: iTyp
    ! Number of atoms
    INTEGER, intent(inout) :: im
    ! Boolean saying if atoms should be kept or not
    LOGICAL, dimension(:), intent(out) :: keep
    ! uAdd(1:3,i : plastic displacement to add to atom i because of Volterra cut
    REAL(kind(0.d0)), dimension(:,:), intent(out) :: uAdd
    ! Periodicity vectors
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    ! Output unit
    INTEGER, intent(in) :: out
    ! Boolean to control verbose mode
    LOGICAL, intent(in) :: verbose

    INTEGER :: nl, il, i, nRemove, nAdd, imNew, imInitial, imm
    LOGICAL :: iCut, vCut
    REAL(kind(0.d0)), dimension(1:3) :: rA, rB, r, xC, b, nAB, uAB, vAB, xpAtom, xp0Atom
    REAL(kind(0.d0)) :: be, z, u, v, rA2, rB2, dotAB, nAB2, factor, Va, vRemove, vAdd
    LOGICAL, dimension(:), allocatable :: added_loop

    INTEGER :: ix, iy, iz
    REAL(kind(0.d0)), dimension(1:3) :: Rcell


    ! Maximal number of atoms allowed by array size
    imm=Size(xp0,2)
    IF ((imm.LT.im).OR.(Size(iTyp).LT.im).OR.(Size(keep).LT.im).OR.(Size(uAdd,2).LT.im)) THEN
            WRITE(0,'(a)') 'Problem with size of arrays on input'
            WRITE(0,'(a,i0)') '  imm = ', imm
            WRITE(0,'(a,i0)') '  im  = ', im 
            STOP '< RemoveAddLoopCut >'
    END IF
    IF ((Size(xp,2).NE.imm).OR.(Size(iTyp).NE.imm).OR.(Size(keep).NE.imm).OR.(Size(uAdd,2).NE.imm)) THEN
            WRITE(0,'(a)') 'Not all arrays have the same size on input'
            STOP '< RemoveAddLoopCut >'
    END IF

    ALLOCATE(added_loop(1:imm))

    ! Initialization
    keep(:)=.TRUE.
    uAdd(:,:)=0.d0

    ! Store initial number of atom
    imNew=im
    imInitial=im

    next_loop: DO nl=1, nLoop

       IF (Debug) WRITE(6,'(a,i0)') "Volterra cut / insertion for loop ", nl
       
       ! Loop center and Burger vector
       xc(:) = cLoop(:,nl)
       b(:) = bLoop(:,nl)

       nRemove=0
       nAdd=0
       vRemove=0.d0
       vAdd=0.d0
       added_loop(:)=.FALSE.    ! tables to record atoms which has been added by this loop

       inside_loop: DO il=1, iLoop(nl)

         ! Triangular loop ABC
         rA(:) = xLoop(:,il-1,nl) - xC(:)
         rB(:) = xLoop(:,il,nl) - xC(:)
         rA2 = Dot_Product(rA, rA)
         rB2 = Dot_Product(rB, rB)
         
         ! Normal to the triangular loop
         nAB = CrossProduct(rA,rB)
         nAB2 = Dot_Product(nAB, nAB)
         IF (nAB2.GT.zero_angle2*rA2*rB2) THEN
                 nAB = nAB/sqrt(nAB2)
         ELSE
                 nAB(:) = 0.d0
         END IF

         ! Edge component of the Burgers vector
         bE = Dot_Product(b,nAB)
         
         IF (Debug) THEN
                 WRITE(6,'(2(a,i0))') "Volterra cut / insertion for portion ", il, " of loop ", nl
                 WRITE(6,'(a,3g13.6)') "  loop corners: A = ", xLoop(:,il-1,nl)
                 WRITE(6,'(a,3g13.6)') "                B = ", xLoop(:,il,nl)
                 WRITE(6,'(a,3g13.6)') "                C = ", xC(:)
                 WRITE(6,'(a,3g13.6)') "       normal:  n = ", nAB(:)
                 WRITE(6,'(a,3g13.6)') "       Burgers: b = ", b(:)
                 WRITE(6,'(a,1g13.6)') "        eddge: be = ", be
         END IF
       
         ! Check that the loop is of V type
         iCut=.FALSE.
         vCut=.FALSE.
         IF (bE.LT.0.d0) THEN
                 ! <0: prismatic loop of interstitial type
                 iCut=.TRUE.
                 ! Volume that is added
                 vAdd = vAdd - 0.5d0*ScalarTripleProduct(rA,rB,b)
                 IF (Debug) WRITE(6,'(a)') "  interstitial loop"
         ELSE IF (bE.gT.0.d0) THEN
                 ! <0: prismatic loop of vacancy type
                 vCut=.TRUE.
                 ! Volume that is removed
                 vRemove = vRemove + 0.5d0*ScalarTripleProduct(rA,rB,b)
                 IF (Debug) WRITE(6,'(a)') "  vacancy loop"
         ELSE IF (bE.EQ.0.d0) THEN
                 ! 0: glissile loop => nothing to be done
                 IF (Debug) WRITE(6,'(a)') "  glissile loop"
                 CYCLE
         END IF


         ! Build projection vectors uAB and vAB:
         !   rA.uAB=1 rB.uAB=0
         !   rA.vAB=0 rB.vAB=1
         dotAB = Dot_Product(rA, rB)

         uAB(:) = dotAB*rB(:) - rB2*rA(:)
         vAB(:) = dotAB*rA(:) - rA2*rB(:)

         factor = 1.d0/( dotAB**2 - rA2*rB2 )
         uAB(:) = factor*uAB(:)
         vAB(:) = factor*vAB(:)

         ! Loop on all image cells
         images_loop: DO ix=-nxImages, nxImages
           DO iy=-nyImages, nyImages
              DO iz=-nzImages, nzImages
                 Rcell(:) = dble(ix)*at(1:3,1) + dble(iy)*at(1:3,2) + dble(iz)*at(1:3,3)
                 IF (DebugPerio) THEN
                         WRITE(6,'(3(a,i0))')  "  periodic cell: ix=", ix, " iy=", iy, " iz=", iz
                         WRITE(6,'(a,3g13.6)') "    Rcell = ", Rcell(1:3)
                 END IF
                 atom_loop: DO i=1, imNew     ! Loop on atoms
                    xp0Atom(:) = xp0(:,i) + Rcell(:)
                    xpAtom(:)  = xp(:,i)  + Rcell(:)
                    r(:) = xp0Atom(:) - xc(:)
                    z = Dot_Product( r(:), nAB(:) )
                    IF ( (z.GT.0.5d0*abs(bE)).OR.(z.LE.-0.5d0*abs(bE)) ) Cycle
                    IF (Debug) WRITE(6,'(a,i0,a)') "  atom ", i, " inside Volterra z-slab"
                    u = Dot_Product( r(:), uAB(:) )
                    IF ( (u.LT.0.d0).OR.(u.GE.1.d0) ) Cycle
                    v = Dot_Product( r(:), vAB(:) )
                    IF ( (v.LT.0.d0).OR.(v.GE.1.d0-u) ) Cycle  
                    IF (Debug) WRITE(6,'(a,i0,a)') "  atom ", i, " inside Volterra cut"
                    IF (vCut) THEN
                            ! Vacancy loop
                            IF (keep(i)) THEN
                                    keep(i)=.FALSE.
                                    nRemove = nRemove+1
                                    IF (Debug) WRITE(out,'(a,i0)') "    remove atom ", i
                            ELSE
                                    IF (Debug) WRITE(out,'(a,i0,a)') "    atom ", i, " already removed"
                            END IF
                    ELSE IF (iCut) THEN
                            ! Interstitial loop: add one atome
                            IF (added_loop(i)) CYCLE ! This atom has already been added
                                                     ! by another portion of the Volterra insertion
                            added_loop(i)=.TRUE.
                            nAdd=nAdd+1
                            im=im+1
                            IF (im.GT.imm) THEN
                                    WRITE(0,'(a)') 'Size of atom arrays is too small'
                                    WRITE(0,'(a)') 'You need to increase value of imm in the input file'
                                    STOP '< RemoveAddLoopCut >'
                            END IF
                            ! Create new atom
                            xp0(1:3,im) = xp0Atom(1:3) - Rcell(1:3)
                            xp(1:3,im)  = xpAtom(1:3)  - Rcell(1:3)
                            iTyp(im) = iTyp(i)
                            keep(im) = keep(i)
                            IF (z.GT.0.d0) THEN
                                    uAdd(:,im) = uAdd(:,i) + 1.d0*b(:)
                            ELSE                            
                                    uAdd(:,im) = uAdd(:,i) - 1.d0*b(:)
                            END IF
                            IF (Debug) THEN
                                    WRITE(out,'(2(a,i0))') "    add atom ", i, " as new atom ", im
                                    WRITE(out,'(a,3g20.12)') "      position: xp0(1:3) = ", xp0(1:3,im)
                                    WRITE(out,'(a,3g20.12)') "      added displacement: uAdd(1:3) = ", uAdd(1:3,im)
                            END IF
                              
                    END IF
                 END DO atom_loop
         END DO ; END DO ; END DO images_loop


       END DO inside_loop

       ! New number of atoms for next loop
       imNew = im

       IF (verbose) THEN
               WRITE(out,'(a,i0,a)') 'Loop ', nl, ' of vacancy / interstitial type'
               WRITE(out,'(2x,2(i0,a))') nRemove, ' atoms removed  -  ',  nAdd, ' atoms added'
               ! Corresponding number of atoms that have been theoretically removed
               Va = Abs( MatDet( at(:,:) ) )/dble(imInitial)
               WRITE(out,'(a,3(f0.2,a))') '  Volume of the Volterra cut, V = ', vRemove, &
                        ', corresponds to ', vRemove/VA , ' atomic volumes ( Va = ', Va,' )'
               WRITE(out,'(a,3(f0.2,a))') '  Volume of the Volterra insertion, V = ', vAdd, &
                        ', corresponds to ', vAdd/VA , ' atomic volumes ( Va = ', Va,' )'
               WRITE(out,*)
       END IF

    END DO next_loop

    DEALLOCATE(added_loop)

       
    IF ( (verbose).AND.(nLoop.GE.2) ) THEN
            WRITE(out,'(a,i0,a)') 'Loop population with ', nLoop, ' loops'
            WRITE(out,'(2x,2(i0,a))') Count(.NOT.keep(1:im)), ' atoms removed  -  ', &
                    im-imInitial, ' atoms added (total)'
            WRITE(out,*)
    END IF

     ! Remove atoms which have been inserted multiple times by dislocation branch cuts
     CALL RemoveDuplicatedAtoms(xp, im, keep, uAdd, out, verbose) 


  END SUBROUTINE RemoveAddLoopCut

  FUNCTION Branch_CutDislo_Cross( x0, x, n) RESULT(cross)
     ! Cross=0, if x0(:) and x(:) are in the same side of the branch cut for dislo n
     ! Cross=+/-1 if x0(:) and x(:) are apart of the branch cut 

     USE disloc_elasticity_ani
     IMPLICIT NONE
     REAL(kind(0.d0)), dimension(1:3), intent(in) :: x0, x
     INTEGER, intent(in) :: n
     REAL(kind(0.d0)) :: cross

     REAL(kind(0.d0)), dimension(1:3) :: s0, s

     ! Default value
     cross=0.d0

     ! Coordinates in the reference frame of the dislocation
     s0(:) = MatMul( rotDislo(:,:,n), x0(:)-cDislo(:,n) )
     s(:) = MatMul( rotDislo(:,:,n), x(:)-cDislo(:,n) )

     IF ((s(1).LE.0.d0).AND.(s0(1).LE.0.d0)) THEN
             IF ( (s0(2).LE.0.d0).AND.(s(2).GT.0.d0) ) THEN
                     cross=1.d0
             ELSE IF ( (s0(2).GT.0.d0).AND.(s(2).LE.0.d0) ) THEN
                     cross=-1.d0
             END IF
     END IF

  END FUNCTION Branch_CutDislo_Cross

  FUNCTION Branch_CutDDipole_Cross( x0, x, n) RESULT(cross)
     ! Cross=0, if x0(:) and x(:) are in the same side of the branch cut for dislo dipole n
     ! Cross=+/-1 if x0(:) and x(:) are apart of the branch cut 

     USE DDipoleModule
     IMPLICIT NONE
     REAL(kind(0.d0)), dimension(1:3), intent(in) :: x0, x
     INTEGER, intent(in) :: n
     REAL(kind(0.d0)) :: cross

     REAL(kind(0.d0)), dimension(1:3) :: s0, s

     ! Default value
     cross=0.d0

     ! Coordinates in the reference frame of the dislocation
     s0(:) = MatMul( rotDDipole(:,:,n), x0(:)-c1DDipole(:,n) )
     s(:) = MatMul( rotDDipole(:,:,n), x(:)-c1DDipole(:,n) )

     IF ((s(1).LE.0.d0).AND.(s0(1).LE.0.d0)) THEN
             IF ( (s0(2).LE.0.d0).AND.(s(2).GT.0.d0) ) THEN
                     cross=1.d0
             ELSE IF ( (s0(2).GT.0.d0).AND.(s(2).LE.0.d0) ) THEN
                     cross=-1.d0
             END IF
     END IF

  END FUNCTION Branch_CutDDipole_Cross

  SUBROUTINE RemoveDuplicatedAtoms(xp, im, keep, uAdd, out, verbose)
     ! Remove atoms which have been duplicated multiple times by dislo branch cuts
     !   variable keep(i) is set to .false. if atom i needs to be removed

     USE Math
     USE Babel_Data, ONLY :  debug, distance_zero2
     IMPLICIT NONE

     ! Atoms position: xp(1:3,i) for atom i (xp0: initial, xp: strained)
     REAL(kind(0.d0)), dimension(:,:), intent(in) :: xp
     ! Number of atoms
     INTEGER, intent(in) :: im
     ! Boolean saying if atoms should be kept or not
     LOGICAL, dimension(:), intent(inout) :: keep
     ! uAdd(1:3,i : plastic displacement to add to atom i because of Volterra cut
     REAL(kind(0.d0)), dimension(:,:), intent(out) :: uAdd
     ! Output unit
     INTEGER, intent(in) :: out
     ! Boolean to control verbose mode
     LOGICAL, intent(in) :: verbose

     INTEGER :: i, j, nKeep0, nKeep1
     REAL(kind(0.d0)) :: d2, du2
     REAL(kind(0.d0)), dimension(1:3) :: dxp, du

     nKeep0=Count(Keep(1:im))

     DO i=1, im
        IF (.NOT.keep(i)) Cycle
        DO j=i+1, im
           IF (.NOT.keep(j)) Cycle
           dxp(:) = xp(:,j) - xp(:,i)
           d2 = SUM( dxp(:)**2 )
           IF (d2.GT.distance_zero2) Cycle
           du(:) = uAdd(:,j) - uAdd(:,i)
           du2 = SUM( du(:)**2 )
           IF (du2.GT.distance_zero2) Cycle
           ! Atoms i and j are duplicated => delete atom j
           keep(j)=.FALSE.
           IF (DEBUG) THEN
                   WRITE(out,'(2(a,i0))') "Duplicated atoms found: ", i, " and ", j
                   WRITE(out,'(a,i0,2(a,3g20.12))') "  atom ", i, ":  position: xp(1:3) = ", xp(1:3,i), &
                        "  added displacement: uAdd(1:3) = ", uAdd(1:3,i)
                   WRITE(out,'(a,i0,2(a,3g20.12))') "  atom ", j, ":  position: xp(1:3) = ", xp(1:3,j), &
                        "  added displacement: uAdd(1:3) = ", uAdd(1:3,j)
           END IF
        END DO
     END DO

     nKeep1=Count(Keep(1:im))

     IF (verbose) THEN
             WRITE(out,*)
             WRITE(out,'(a,i0)') 'Volterra cut: number of atoms which have been deleted because of multiple insertion: ', &
                     nKeep0-nKeep1
             WRITE(out,*)
     END IF

  END SUBROUTINE RemoveDuplicatedAtoms

END MODULE Volterra_Module
