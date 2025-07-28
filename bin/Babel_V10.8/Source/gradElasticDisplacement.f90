MODULE GradElasticDisplacementModule

   IMPLICIT NONE

   SAVE
   
   ! Method used to select pattern for each atom: patternSelectionMethod
   !   = 1 : Select pattern which maximises the number of matchs 
   !         look to the NORM OF THE DISPLACEMENT if 2 patterns are equivalent
   !   = 2 : Select pattern which maximises the number of matchs 
   !         look to the norm of the ELASTIC GRADIENT if 2 patterns are equivalent
   !   = 3 : Select pattern which minimizes the number of unmatched bonds 
   !         look to the NORM OF THE DISPLACEMENT if 2 patterns are equivalent
   !   = 4 : Select pattern which minimizes the number of unmatched bonds 
   !         look to the norm of the ELASTIC GRADIENT if 2 patterns are equivalent
   ! (default: 3)
   INTEGER :: patternSelectionMethod

   ! The pattern is read in the file patternFile
   CHARACTER(len=200) :: patternFile

   ! Verbosity for printing information
   INTEGER, parameter, private :: verbosity_max=4   
   ! Verbosity needed for printing debug information for pattern search
   INTEGER, parameter, private :: verbosity_debug_pattern=20
   ! Verbosity needed for printing debug information for calculation of the
   !   gradient of the elastic displacement
   INTEGER, parameter, private :: verbosity_debug_gradElasticDisplacement=20

   ! Use periodic boundary conditions to calculate elastic displacement gradient ...
   LOGICAL, parameter, private :: perioElasticDisplacement=.true.

   ! Gradient of elastic displacement for each atom
   !  gradElasticDisplacement(i,j, n) = d U_i / d X_j for atom n
   REAL(kind(0.d0)), dimension(:,:,:), allocatable :: gradElasticDisplacement


   ! Number of different atomic patterns in reference structure (and maximal value)
   INTEGER, parameter, private :: max_nPat=2000
   INTEGER, private :: nPat=0

   ! Maximal number of neighbours in a pattern
   INTEGER, private :: max_nNeighPat=0

   ! Atomic patterns definition:
   !   nNeighPat(n): number of neighbours in pattern n
   INTEGER, dimension(1:max_nPat), private :: nNeighPat
   !   pattern(1:3,i,n): coordinate of neighbour i for pattern n
   REAL(kind(0.d0)), dimension(:,:,:), allocatable, private :: pattern
   !   corresponding normalized vectors
   REAL(kind(0.d0)), dimension(:,:,:), allocatable, private :: uPattern
   !   pattern nickname
   CHARACTER(len=50), dimension(1:max_nPat), private :: labelPattern

   ! refPattern(n) and inpPattern(n): pattern number for atom n in reference and input structures
   INTEGER, dimension(:), allocatable :: refPattern
   INTEGER, dimension(:), allocatable :: inpPattern

   ! Tolerance for atomic distances (relative)
   REAL(kind(0.d0)), parameter, private :: tolerance=1.d-1
   REAL(kind(0.d0)), parameter, private :: tolerance2=tolerance**2

   ! Threshold angle for pairing vectors
   REAL(kind(0.d0)) :: patternAngleThreshold=10.d0      ! degrees

   PRIVATE :: mat3inv

   Private :: PatternMatch

CONTAINS

  SUBROUTINE InitRefPattern(im0, xp0, at0, nNeigh, iNeigh, out)
     ! Record atomic patterns found in reference structure

     USE Babel_data, ONLY : verbosity, imm
     IMPLICIT NONE

     ! Number of atoms
     INTEGER, intent(in) :: im0
     ! Atoms coordinates in reference configuration
     REAL(kind(0.d0)), dimension(:,:), intent(in) :: xp0
     ! Periodicity vectors of the reference structure: at0(1:3,1), at0(1:3,2), ...
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at0
     ! Number of neighbours for each atom 
     INTEGER, dimension(:), intent(in) :: nNeigh
     ! Neighbour indexes for each atom i: iNeigh(1:nNeigh(i), i)
     INTEGER, dimension(:,:), intent(in) :: iNeigh
     ! Output unit 
     INTEGER, intent(in), optional :: out

     ! Local variables
     INTEGER :: i, n, j, nMax, np
     LOGICAL :: new_pattern
     REAL(kind(0.d0)) :: invR0
     REAL(kind(0.d0)), dimension(1:3) :: ds0
     REAL(kind(0.d0)), dimension(:,:), allocatable :: R0
     REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at0
     ! Number of atoms in pattern p: nAtoms_in_pat(p)
     INTEGER, dimension(1:max_nPat) :: nAtoms_in_pat

     IF (verbosity.GE.verbosity_max) THEN
             WRITE(out,*)
             WRITE(out,'(a)') 'Look for pattern definition in reference structure'
             WRITE(out,*)
     END IF

     !-----------------------------------
     ! table recording pattern number for atom n in reference structure 
     IF (Allocated(refPattern)) Deallocate(refPattern)
     Allocate(refPattern(1:imm))
     refPattern(:)=0

     ! Maximal number of neighbours
     max_nNeighPat = MaxVal(nNeigh(1:im0))

     ! Allocation and initialization for pattern definition
     nPat=0
     IF (Allocated(pattern)) Deallocate(pattern)
     ALLOCATE(pattern(1:3,1:max_nNeighPat,1:max_nPat))
     IF (Allocated(uPattern)) Deallocate(uPattern)
     ALLOCATE(uPattern(1:3,1:max_nNeighPat,1:max_nPat))
     pattern(:,:,:)=0.d0
     uPattern(:,:,:)=0.d0
     nNeighPat(:)=0
     labelPattern(:)=''

     ! table counting number of atoms in each pattern
     nAtoms_in_pat(:) = 0

     ! Local array
     ALLOCATE(R0(1:3,1:max_nNeighPat))

     ! Inverse periodicity vectors
     IF (perioElasticDisplacement) THEN
             CALL Mat3Inv(at0, inv_at0)
     ELSE
             inv_at0(:,:)=0.d0
     END IF

     !-----------------------------------
     ! Loop on all atoms of reference structure
     atom_loop: DO i=1, im0
        
        ! initialization
        new_pattern=.TRUE.

        ! number of neighbours
        nMax = nNeigh(i)

        ! Loop on neighbours
        neighbour_loop: DO n=1, nMax
          ! Index of neighbour atom
          j = iNeigh(n, i)
          R0(1:3,n) = xp0(1:3,j) - xp0(1:3,i)
          ! Apply periodic boundary conditions to atom coordinates
          ! (works only ! for ~orthogonal vectors)
          IF (perioElasticDisplacement) THEN
                  ds0(:) = MatMul( inv_at0(:,:), R0(1:3,n) )
                  ds0(:) = aNInt(ds0(:))
                  R0(1:3,n) = R0(1:3,n) - MatMul( at0(:,:), ds0(:) )
          END IF
        END DO neighbour_loop

        IF (verbosity.GE.verbosity_debug_pattern) THEN
                WRITE(out,*)
                WRITE(out,'(a)') 'Potential new pattern definition'
                WRITE(out,'(a,i0)') '  number of neighbours in pattern definition: ', nMax
                DO n=1, nMax
                   WRITE(out,'(a,i3,a,3g14.6)')'    atom ', n, ': R(1:3) = ', R0(1:3,n)
                END DO
        END IF

        ! loop on patterns
        pattern_loop: DO np=1, nPat
           IF (verbosity.GE.verbosity_debug_pattern) &
                   WRITE(out,'(a,i0)') 'Compare potential new pattern with pattern ', np
           IF ( PatternMatch( R0(:,:), nMax, pattern(:,:,np), nNeighPat(np) ) ) THEN
                   new_pattern=.false.  ! Both patterns are equivalent
                   EXIT pattern_loop
           END IF
        END DO pattern_loop

        ! new pattern has been detectect
        IF (new_pattern) THEN
                nPat = nPat + 1
                IF (nPat.GT.max_nPat) THEN
                        WRITE(0,'(a)') 'Too many patterns detected'
                        WRITE(0,'(a)') 'You need to modify value for max_nPat&
                               & in file elasticStrain.f90 and to recompile'
                        STOP '< InitRefPatterns >'
                END IF
                nNeighPat(nPat) = nMax          ! Number of neighbours
                pattern(:,:,nPat) = R0(:,:)     ! Neighbour positions
                ! Normalized neighbour vectors
                DO n=1, nMax
                   invR0 = 1.d0/sqrt( Sum( R0(1:3,n)**2 ) )
                   uPattern(:,n,nPat) = R0(:,n)*invR0
                END DO

                ! Atom i corresponds to pattern nPat
                refPattern(i) = nPat

                ! Pattern np counts one more atom
                nAtoms_in_pat(np) = nAtoms_in_pat(np) + 1

                IF (verbosity.GE.verbosity_debug_pattern) THEN
                        WRITE(out,*)
                        WRITE(out,'(a)') 'New pattern found'
                        CALL PrintPattern(nPat,out)
                        WRITE(out,*)
                END IF
        ELSE 
                ! Atom i corresponds to pattern np
                refPattern(i) = np

                ! Pattern np counts one more atom
                nAtoms_in_pat(np) = nAtoms_in_pat(np) + 1

                IF (verbosity.GE.verbosity_debug_pattern) &
                        WRITE(out,'(a,i0)') 'Pattern already exists: equivalent to pattern ', np
        END IF

     END DO atom_loop

     Deallocate(R0)

     IF (verbosity.GE.verbosity_max) THEN
             CALL PrintAllPatterns(out)
             WRITE(out,*)
             DO np=1, nPat
                WRITE(out,'(2(a,i0))') 'Number of atoms corresponding to pattern ', np,': ', &
                        nAtoms_in_pat(np)
             END DO
             WRITE(out,*)
     END IF
     IF (verbosity.GE.verbosity_debug_pattern) THEN
             WRITE(out,*)
             DO i=1, im0
                WRITE(out,'(2(a,i0))') ' atom ', i, ' :  pattern = ', refPattern(i)
             END DO
     END IF

  END SUBROUTINE InitRefPattern

  FUNCTION PatternMatch( rPat1, n1, rPat2, n2) RESULT(patMatch)
     ! Return true if patterns defined by rPat1(1:3,:) and rPat2(1:3,:) are equivalent

     IMPLICIT NONE

     ! Number of neighbours in each pattern
     INTEGER, intent(in) :: n1, n2
     ! rPat1(1:3,i) coordinates of vector i
     REAL(kind(0.d0)), dimension(:,:), intent(in) :: rPat1, rPat2
     LOGICAL :: patMatch

     ! Local variables
     INTEGER :: i1, i2
     REAL(kind(0.d0)), dimension(1:3) :: R1, R2
     REAL(kind(0.d0)) :: inv_R1square
     LOGICAL :: atomMatch
     LOGICAL, dimension(1:n2) :: usedPat2

     ! Check that both patterns have the same number of neighbours
     IF (n1.NE.n2) THEN
             patMatch=.false.
             RETURN
     ELSE
             patMatch=.true.
     END IF

     ! Check atom positions
     atomPattern1_loop: DO i1=1, n1

        ! initialization
        R1(:) = rPat1(:,i1)
        atomMatch=.false.
        usedPat2(:)=.false.

        inv_R1square = 1.d0/Sum( R1(:)**2 )
     
        atomPattern2_loop: DO i2=1, n1

           ! Check if atom in pattern 2 has not been already used
           IF (usedPat2(i2)) Cycle

           R2(:) = rPat2(:,i2)

           IF ( Sum( (R2(:)-R1(:))**2 )*inv_R1square.LE.tolerance2 ) THEN
                   ! The two atoms match
                   atomMatch=.true.
                   EXIT atomPattern2_loop
           END IF

        END DO atomPattern2_loop

        IF (atomMatch) THEN
                ! rPat1(:,i1) and rPat(:,i2) are equivalent
                usedPat2(i2)=.true.
        ELSE
                ! did not manage to find an atom corresponding to rPat1(:,i1) in pattern 2
                patMatch=.false.
                RETURN
        END IF

     END DO atomPattern1_loop

  END FUNCTION PatternMatch

  SUBROUTINE PrintAllPatterns(out)
    ! Print pattern definition on output unit out

    IMPLICIT NONE

    INTEGER, intent(in) :: out
    INTEGER :: n
    
    WRITE(out,*)
    WRITE(out,'(a,i0)') 'Number of patterns, nPat = ', nPat
    DO n=1, nPat
       CALL PrintPattern(n, out)
    END DO
    WRITE(out,*)

  END SUBROUTINE PrintAllPatterns

  SUBROUTINE PrintPattern(n, out)
    ! Print pattern definition on output unit out for pattern n

    IMPLICIT NONE

    INTEGER, intent(in) :: n
    INTEGER, intent(in) :: out
    INTEGER :: i
    
    WRITE(out,*)
    WRITE(out,'(a,i0)') 'Pattern ', n
    IF (Trim(labelPattern(n)).NE.'') WRITE(out,'(a)') Trim(labelPattern(n))
    WRITE(out,'(a,i0)') '  number of neighbours in pattern definition: ', nNeighPat(n)
    DO i=1, nNeighPat(n)
       WRITE(out,'(a,i3,a,3g14.6)')'    atom ', i, ': R(1:3) = ', pattern(1:3,i,n)
    END DO

  END SUBROUTINE PrintPattern

  SUBROUTINE RotatePattern(rot)
    ! Rotate all patterns with rotation matrix rot

    IMPLICIT NONE

    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_rot, temp
    INTEGER :: i, n

   ! Identity matrix
    REAL(kind(0.d0)), dimension(1:3,1:3), parameter :: matId=ReShape( &
        (/ 1.d0, 0.d0, 0.d0, &
           0.d0, 1.d0, 0.d0, &
           0.d0, 0.d0, 1.d0 /), (/ 3,3 /) )

    ! Check that rot is a rotation matrix
    inv_rot = Transpose( rot )
    temp = MatMul( inv_rot, rot ) - matId
    IF ( Sum( temp(:,:)**2 ).GT.1.d-12) THEN
            WRITE(0,'(a)') 'Rot(1:3,1:3) is not a unitary matrix'
            STOP '< RotatePattern >'
    END IF

    DO n=1, nPat
       DO i=1, nNeighPat(n)
           pattern(:,i,n) = MatMul( rot(:,:),  pattern(:,i,n) )
          uPattern(:,i,n) = MatMul( rot(:,:), uPattern(:,i,n) )
       END DO
    END DO

  END SUBROUTINE RotatePattern

  SUBROUTINE ReadPattern(inp, out)
     ! Read pattern definition in file connected to unit inp

     USE Babel_data, ONLY : verbosity
     IMPLICIT NONE
     ! input and output units
     INTEGER, intent(in) :: inp, out

     INTEGER :: i, n, nNeigh

     ! Number of different atomic patterns in reference structure
     CALL Comment(inp)
     READ(inp,*) nPat
     IF ( (nPat.LE.0) .OR. (nPat.GT.max_nPat) ) THEN
             WRITE(0,'(a,i0)') 'nPat = ', nPat
             STOP '< ReadPattern >'
     END IF

     ! First pass to determine the size needed for the arrays
     max_nNeighPat=0
     DO n=1, nPat
        ! number of neighbours in pattern n
        CALL Comment(inp)
        READ(inp,*) nNeigh
        IF (nNeigh.GT.max_nNeighPat) max_nNeighPat=nNeigh
        ! coordinate of neighbour i for pattern n
        DO i=1, nNeigh
           CALL Comment(inp)
           READ(inp,*) 
        END DO
     END DO
        
     REWIND(inp)

     ! Allocation and initialization for pattern definition
     IF (Allocated(pattern)) Deallocate(pattern)
     ALLOCATE(pattern(1:3,1:max_nNeighPat,1:max_nPat))
     IF (Allocated(uPattern)) Deallocate(uPattern)
     ALLOCATE(uPattern(1:3,1:max_nNeighPat,1:max_nPat))
     pattern(:,:,:)=0.d0
     uPattern(:,:,:)=0.d0
     nNeighPat(:)=0
     labelPattern(:)=''

     ! Number of different atomic patterns in reference structure
     CALL Comment(inp)
     READ(inp,*) nPat

     DO n=1, nPat
        ! nNeighPat(n): number of neighbours in pattern n
        CALL Comment(inp)
        READ(inp,*) nNeighPat(n)
        ! pattern(1:3,i,n): coordinate of neighbour i for pattern n
        DO i=1, nNeighPat(n)
           CALL Comment(inp)
           READ(inp,*) pattern(1:3,i,n)
           uPattern(1:3,i,n) = pattern(1:3,i,n) / sqrt( Sum( pattern(1:3,i,n)**2) )
        END DO
     END DO

     IF (verbosity.GE.verbosity_max) THEN
             WRITE(out,*)
             WRITE(out,'(a)') 'Pattern read in file'
             CALL PrintAllPatterns(out)
     END IF

  END SUBROUTINE ReadPattern

  SUBROUTINE WritePattern(out)
     ! Write pattern definition in file connected to unit inp
     ! The file can be then read with the subroutine ReadPattern

     IMPLICIT NONE
     ! Output unit
     INTEGER, intent(in) :: out

     INTEGER :: n, i

     WRITE(out,'(a)') "# Number of different patterns"
     WRITE(out,'(i0)') nPat
     WRITE(out,*)

     DO n=1, nPat
        WRITE(out,'(a,i0)') "# Definition of pattern ", n
        IF (Trim(labelPattern(n)).NE.'') WRITE(out,'(2a)') "# ", Trim(labelPattern(n))
        WRITE(out,'(a)')    "#    number of neighbours in pattern"
        WRITE(out,'(i0)')   nNeighPat(n)
        WRITE(out,'(a)')    "#    coordinates of the different neighbours"
        DO i=1, nNeighPat(n)
           WRITE(out,'(3g20.12)') pattern(1:3,i,n)
        END DO
        WRITE(out,*)
     END DO

  END SUBROUTINE WritePattern

  SUBROUTINE InitHcpPattern(alat,coa)
     ! Create pattern corresponding to hcp structure

     IMPLICIT NONE

     ! Lattice parameters and c/a ratio of hcp lattice
     REAL(kind(0.d0)), intent(in) :: alat, coa

     INTEGER :: i, n

     ! Allocation and initialization for pattern definition
     max_nNeighPat = 13
     IF (Allocated(pattern)) Deallocate(pattern)
     ALLOCATE(pattern(1:3,1:max_nNeighPat,1:max_nPat))
     IF (Allocated(uPattern)) Deallocate(uPattern)
     ALLOCATE(uPattern(1:3,1:max_nNeighPat,1:max_nPat))
     pattern(:,:,:)=0.d0
     uPattern(:,:,:)=0.d0
     nNeighPat(:)=0
     labelPattern(:) = ''
     
     ! Number of different atomic patterns in reference structure
     nPat=2

     ! Pattern 1:
     !   nickname
     labelPattern(1) = "hcp - plane A"
     !   number of neighbours
     nNeighPat(1) = 12
     !   coordinates of each neighbour
     pattern(1:3, 1,1) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 2,1) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 3,1) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 4,1) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 5,1) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 6,1) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 7,1) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3, 8,1) = alat*(/   0.d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
     pattern(1:3, 9,1) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3,10,1) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)
     pattern(1:3,11,1) = alat*(/   0.d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern(1:3,12,1) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)

     ! Pattern 2:
     !   nickname
     labelPattern(2) = "hcp - plane B"
     !   number of neighbours
     nNeighPat(2) = 12
     !   coordinates of each neighbour
     pattern(1:3, 1,2) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 2,2) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 3,2) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 4,2) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 5,2) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 6,2) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 7,2) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3, 8,2) = alat*(/   0.d0,  sqrt(3.d0)/3.d0,  0.5d0*coa /)
     pattern(1:3, 9,2) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3,10,2) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)
     pattern(1:3,11,2) = alat*(/   0.d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern(1:3,12,2) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)

     ! Reduced coordinates for each pattern
     DO n=1, nPat
        DO i=1, nNeighPat(n)
           uPattern(1:3,i,n) = pattern(1:3,i,n) / sqrt( Sum( pattern(1:3,i,n)**2) )
        END DO
     END DO

  END SUBROUTINE InitHcpPattern

  SUBROUTINE InitBccPattern(alat)
     ! Create pattern corresponding to bcc structure

     IMPLICIT NONE

     ! Lattice parameters of bcc lattice
     REAL(kind(0.d0)), intent(in) :: alat

     INTEGER :: i, n

     ! Allocation and initialization for pattern definition
     max_nNeighPat = 8
     IF (Allocated(pattern)) Deallocate(pattern)
     ALLOCATE(pattern(1:3,1:max_nNeighPat,1:max_nPat))
     IF (Allocated(uPattern)) Deallocate(uPattern)
     ALLOCATE(uPattern(1:3,1:max_nNeighPat,1:max_nPat))
     pattern(:,:,:)=0.d0
     uPattern(:,:,:)=0.d0
     nNeighPat(:)=0
     labelPattern(:) = ''
     
     ! Number of different atomic patterns in reference structure
     nPat=1

     ! Pattern 1:
     !   nickname
     labelPattern(1) = "bcc"
     !   number of neighbours
     nNeighPat(1) = 8
     !   coordinates of each neighbour
     pattern(1:3, 1,1) = alat*(/ -sqrt(6.d0)/3.d0,   0.d0          ,  sqrt(3.d0)/6.d0 /)
     pattern(1:3, 2,1) = alat*(/  sqrt(6.d0)/6.d0,  sqrt(2.d0)/2.d0,  sqrt(3.d0)/6.d0 /)
     pattern(1:3, 3,1) = alat*(/  sqrt(6.d0)/6.d0, -sqrt(2.d0)/2.d0,  sqrt(3.d0)/6.d0 /)
     pattern(1:3, 4,1) = alat*(/   0.d0          ,   0.d0          ,  sqrt(3.d0)      /)
     pattern(1:3, 5,1) = alat*(/   0.d0          ,   0.d0          , -sqrt(3.d0)      /)
     pattern(1:3, 6,1) = alat*(/  sqrt(6.d0)/3.d0,   0.d0          , -sqrt(3.d0)/6.d0 /)
     pattern(1:3, 7,1) = alat*(/ -sqrt(6.d0)/6.d0,  sqrt(2.d0)/2.d0, -sqrt(3.d0)/6.d0 /)
     pattern(1:3, 8,1) = alat*(/ -sqrt(6.d0)/6.d0, -sqrt(2.d0)/2.d0, -sqrt(3.d0)/6.d0 /)

     ! Reduced coordinates for each pattern
     DO n=1, nPat
        DO i=1, nNeighPat(n)
           uPattern(1:3,i,n) = pattern(1:3,i,n) / sqrt( Sum( pattern(1:3,i,n)**2) )
        END DO
     END DO

  END SUBROUTINE InitBccPattern

  SUBROUTINE HcpPattern_AddBasalFault(alat,coa)
     ! Create pattern corresponding to basal stacking faults in hcp structure

     IMPLICIT NONE

     ! Lattice parameters and c/a ratio of hcp lattice
     REAL(kind(0.d0)), intent(in) :: alat, coa

     INTEGER :: i, n

     ! Allocation and initialization for pattern definition
     IF ( ( .NOT. Allocated(pattern) ).OR.( .NOT. Allocated(uPattern) ).OR.( nPat.LE.0 ) ) THEN
             WRITE(0,'(a)') 'You need to initialize HCP pattern before'
             STOP '< HcpPattern_AddBasalFault >'
     END IF

     IF (max_nNeighPat.LT.12) THEN
             WRITE(0,'(a,i0)') 'Problem with maximal number of neighbours in a pattern :  max_nNeighPat = ', max_nNeighPat
             WRITE(0,'(a)')    'It should be at least equal to 12'
             STOP '< HcpPattern_AddBasalFault >'
     END IF
     
     ! Pattern 1:
     !   nickname
     labelPattern(nPat+1) = "fcc - stacking ABC"
     !   number of neighbours
     nNeighPat(nPat+1) = 12
     !   coordinates of each neighbour
     pattern(1:3, 1,nPat+1) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 2,nPat+1) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 3,nPat+1) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 4,nPat+1) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 5,nPat+1) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 6,nPat+1) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 7,nPat+1) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3, 8,nPat+1) = alat*(/   0.d0,  sqrt(3.d0)/3.d0,  0.5d0*coa /)
     pattern(1:3, 9,nPat+1) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3,10,nPat+1) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)
     pattern(1:3,11,nPat+1) = alat*(/   0.d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern(1:3,12,nPat+1) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)

     ! Pattern 2:
     !   nickname
     labelPattern(nPat+2) = "fcc - stacking ACB"
     !   number of neighbours
     nNeighPat(nPat+2) = 12
     !   coordinates of each neighbour
     pattern(1:3, 1,nPat+2) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 2,nPat+2) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 3,nPat+2) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 4,nPat+2) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 5,nPat+2) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 6,nPat+2) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 7,nPat+2) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3, 8,nPat+2) = alat*(/   0.d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
     pattern(1:3, 9,nPat+2) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3,10,nPat+2) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)
     pattern(1:3,11,nPat+2) = alat*(/   0.d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern(1:3,12,nPat+2) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)

     ! Reduced coordinates for each pattern
     DO n=nPat+1, nPat+2
        DO i=1, nNeighPat(n)
           uPattern(1:3,i,n) = pattern(1:3,i,n) / sqrt( Sum( pattern(1:3,i,n)**2) )
        END DO
     END DO

     ! Number of different atomic patterns in reference structure
     nPat = nPat+2

  END SUBROUTINE HcpPattern_AddBasalFault

  SUBROUTINE HcpPattern_AddPrismFault(alat,coa)
     ! Create pattern corresponding to prismatic stacking faults in hcp structure
     !   fault plane: (10-10)
     !   fault vector: 1/6 [1-210]

     IMPLICIT NONE

     ! Lattice parameters and c/a ratio of hcp lattice
     REAL(kind(0.d0)), intent(in) :: alat, coa

     INTEGER :: i, n

     ! Allocation and initialization for pattern definition
     IF ( ( .NOT. Allocated(pattern) ).OR.( .NOT. Allocated(uPattern) ).OR.( nPat.LE.0 ) ) THEN
             WRITE(0,'(a)') 'You need to initialize HCP pattern before'
             STOP '< HcpPattern_AddPrismFault >'
     END IF

     IF (max_nNeighPat.LT.13) THEN
             WRITE(0,'(a,i0)') 'Problem with maximal number of neighbours in a pattern :  max_nNeighPat = ', max_nNeighPat
             WRITE(0,'(a)')    'It should be at least equal to 13'
             STOP '< HcpPattern_AddPrismFault >'
     END IF
     
     ! Pattern 1:
     !   nickname
     labelPattern(nPat+1) = "hcp - prism fault A"
     !   number of neighbours
     nNeighPat(nPat+1) = 13
     !   coordinates of each neighbour
     pattern(1:3, 1,nPat+1) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 2,nPat+1) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 3,nPat+1) = alat*(/  0.d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 4,nPat+1) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 5,nPat+1) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 6,nPat+1) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3, 7,nPat+1) = alat*(/  0.5d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
     pattern(1:3, 8,nPat+1) = alat*(/ -0.5d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
     pattern(1:3, 9,nPat+1) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3,10,nPat+1) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)
     pattern(1:3,11,nPat+1) = alat*(/  0.5d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern(1:3,12,nPat+1) = alat*(/ -0.5d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern(1:3,13,nPat+1) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)

     ! Pattern 2:
     !   nickname
     labelPattern(nPat+2) = "hcp - prism fault B"
     !   number of neighbours
     nNeighPat(nPat+2) = 13
     !   coordinates of each neighbour
     pattern(1:3, 1,nPat+2) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 2,nPat+2) = alat*(/   0.d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 3,nPat+2) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 4,nPat+2) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern(1:3, 5,nPat+2) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
     pattern(1:3, 6,nPat+2) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3, 7,nPat+2) = alat*(/   0.5d0,  sqrt(3.d0)/3.d0, 0.5d0*coa /)
     pattern(1:3, 8,nPat+2) = alat*(/  -0.5d0,  sqrt(3.d0)/3.d0, 0.5d0*coa /)
     pattern(1:3, 9,nPat+2) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern(1:3,10,nPat+2) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)
     pattern(1:3,11,nPat+2) = alat*(/  0.5d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern(1:3,12,nPat+2) = alat*(/ -0.5d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern(1:3,13,nPat+2) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)


     ! Reduced coordinates for each pattern
     DO n=nPat+1, nPat+2
        DO i=1, nNeighPat(n)
           uPattern(1:3,i,n) = pattern(1:3,i,n) / sqrt( Sum( pattern(1:3,i,n)**2) )
        END DO
     END DO

     ! Number of different atomic patterns in reference structure
     nPat = nPat+2

  END SUBROUTINE HcpPattern_AddPrismFault

  SUBROUTINE HcpPattern_AddPy1Fault(alat,coa)
     ! Create pattern corresponding to prismatic stacking faults in hcp structure
     !   fault plane: (10-10)
     !   fault vector: 1/6 [1-210]

     IMPLICIT NONE

     ! Lattice parameters and c/a ratio of hcp lattice
     REAL(kind(0.d0)), intent(in) :: alat, coa

     INTEGER :: i, n

     ! Allocation and initialization for pattern definition
     IF ( ( .NOT. Allocated(pattern) ).OR.( .NOT. Allocated(uPattern) ).OR.( nPat.LE.0 ) ) THEN
             WRITE(0,'(a)') 'You need to initialize HCP pattern before'
             STOP '< HcpPattern_AddPy1Fault >'
     END IF

     IF (max_nNeighPat.LT.12) THEN
             WRITE(0,'(a,i0)') 'Problem with maximal number of neighbours in a pattern :  max_nNeighPat = ', max_nNeighPat
             WRITE(0,'(a)')    'It should be at least equal to 12'
             STOP '< HcpPattern_AddPy1Fault >'
     END IF
     
     !=========
     ! Full solution
     !!$ CALL HcpPattern_AddPy1Twin(alat, coa, 1)
     !!$ !CALL HcpPattern_AddPy1Twin(alat, coa, 5)
     !!$ labelPattern(nPat-1) = "hcp - (-1011) pyramidal fault, plane A"
     !!$ labelPattern(nPat)   = "hcp - (-1011) pyramidal fault, plane B"
     !!$ CALL HcpPattern_AddPy1Twin(alat, coa, 7)
     !!$ !CALL HcpPattern_AddPy1Twin(alat, coa, 11)
     !!$ labelPattern(nPat-1) = "hcp - (10-11) pyramidal fault, plane A"
     !!$ labelPattern(nPat)   = "hcp - (10-11) pyramidal fault, plane B"

     !=========
     ! Solution which only consider displacement along the [1-210] direction
     !
      ! Pattern 1:
      !   nickname
      labelPattern(nPat+1) = "hcp - (-1011) pyramidal fault A"
      !   number of neighbours
      nNeighPat(nPat+1) = 12
      !   coordinates of each neighbour
      pattern(1:3, 1,nPat+1) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
      pattern(1:3, 2,nPat+1) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 3,nPat+1) = alat*(/   0.d0, -sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 4,nPat+1) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 5,nPat+1) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
      pattern(1:3, 6,nPat+1) = alat*(/   0.d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
      pattern(1:3, 7,nPat+1) = alat*(/  0.5d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
      pattern(1:3, 8,nPat+1) = alat*(/ -0.5d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
      pattern(1:3, 9,nPat+1) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)
      pattern(1:3,10,nPat+1) = alat*(/  0.5d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)
      pattern(1:3,11,nPat+1) = alat*(/ -0.5d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)
      pattern(1:3,12,nPat+1) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)

      ! Pattern 2:
      !   nickname
      labelPattern(nPat+2) = "hcp - (-1011) pyramidal fault B"
      !   number of neighbours
      nNeighPat(nPat+2) = 12
      !   coordinates of each neighbour
      pattern(1:3, 1,nPat+2) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
      pattern(1:3, 2,nPat+2) = alat*(/   0.d0,  sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 3,nPat+2) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 4,nPat+2) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 5,nPat+2) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
      pattern(1:3, 6,nPat+2) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
      pattern(1:3, 7,nPat+2) = alat*(/  0.5d0,  sqrt(3.d0)/3.d0,  0.5d0*coa /)
      pattern(1:3, 8,nPat+2) = alat*(/ -0.5d0,  sqrt(3.d0)/3.d0,  0.5d0*coa /)
      pattern(1:3, 9,nPat+2) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
      pattern(1:3,10,nPat+2) = alat*(/   0.d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)
      pattern(1:3,11,nPat+2) = alat*(/  0.5d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)
      pattern(1:3,12,nPat+2) = alat*(/ -0.5d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)

      ! Pattern 3:
      !   nickname
      labelPattern(nPat+3) = "hcp - (10-11) pyramidal fault A"
      !   number of neighbours
      nNeighPat(nPat+3) = 12
      !   coordinates of each neighbour
      pattern(1:3, 1,nPat+3) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
      pattern(1:3, 2,nPat+3) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 3,nPat+3) = alat*(/   0.d0, -sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 4,nPat+3) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 5,nPat+3) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
      pattern(1:3, 6,nPat+3) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
      pattern(1:3, 7,nPat+3) = alat*(/  0.5d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
      pattern(1:3, 8,nPat+3) = alat*(/ -0.5d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
      pattern(1:3, 9,nPat+3) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
      pattern(1:3,10,nPat+3) = alat*(/   0.d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)
      pattern(1:3,11,nPat+3) = alat*(/  0.5d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)
      pattern(1:3,12,nPat+3) = alat*(/ -0.5d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)

      ! Pattern 4:
      !   nickname
      labelPattern(nPat+4) = "hcp - (10-11) pyramidal fault B"
      !   number of neighbours
      nNeighPat(nPat+4) = 12
      !   coordinates of each neighbour
      pattern(1:3, 1,nPat+4) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
      pattern(1:3, 2,nPat+4) = alat*(/   0.d0,  sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 3,nPat+4) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 4,nPat+4) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
      pattern(1:3, 5,nPat+4) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
      pattern(1:3, 6,nPat+4) = alat*(/   0.d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
      pattern(1:3, 7,nPat+4) = alat*(/  0.5d0,  sqrt(3.d0)/3.d0,  0.5d0*coa /)
      pattern(1:3, 8,nPat+4) = alat*(/ -0.5d0,  sqrt(3.d0)/3.d0,  0.5d0*coa /)
      pattern(1:3, 9,nPat+4) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)
      pattern(1:3,10,nPat+4) = alat*(/  0.5d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)
      pattern(1:3,11,nPat+4) = alat*(/ -0.5d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)
      pattern(1:3,12,nPat+4) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)

      ! Reduced coordinates for each pattern
      DO n=nPat+1, nPat+4
         DO i=1, nNeighPat(n)
            uPattern(1:3,i,n) = pattern(1:3,i,n) / sqrt( Sum( pattern(1:3,i,n)**2) )
         END DO
      END DO

      ! Number of different atomic patterns in reference structure
      nPat = nPat+4

  END SUBROUTINE HcpPattern_AddPy1Fault

  SUBROUTINE HcpPattern_AddPy1Twin(alat,coa, nSys)
     ! Create pattern corresponding to 1st order pyramidal twins of hcp structure
     !   if nSys is given, only the pattern corresponding to the system nSys is !   added

     IMPLICIT NONE

     ! Lattice parameters and c/a ratio of hcp lattice
     REAL(kind(0.d0)), intent(in) :: alat, coa
     ! System index (optional)
     INTEGER, intent(in), optional :: nSys

     INTEGER :: i, n, n1, n2, nSys1, nSys2
     INTEGER, dimension(1:4,1:12) :: nPy
     REAL(kind(0.d0)), dimension(4,3) :: aHCP
     REAL(kind(0.d0)) :: inv_norm
     REAL(kind(0.d0)), dimension(1:3) :: normal
     REAL(kind(0.d0)), dimension(1:3,1:12,1:2) :: pattern0

     ! Allocation and initialization for pattern definition
     IF ( ( .NOT. Allocated(pattern) ).OR.( .NOT. Allocated(uPattern) ).OR.( nPat.LE.0 ) ) THEN
             WRITE(0,'(a)') 'You need to initialize HCP pattern before'
             STOP '< HcpPattern_AddPy1Twin >'
     END IF

     IF (max_nNeighPat.LT.12) THEN
             WRITE(0,'(a,i0)') 'Problem with maximal number of neighbours in a pattern :  max_nNeighPat = ', max_nNeighPat
             WRITE(0,'(a)')    'It should be at least equal to 12'
             STOP '< HcpPattern_AddPy1Twin >'
     END IF
     
     ! Cartesian coordinates of the hcp lattice
     !   with factor for normal to planes.
     aHCP(1,1:3) = (/   1.d0,       0.d0,       0.d0 /)
     aHCP(2,1:3) = (/ -0.5d0,  sqrt(3.d0)/2.d0, 0.d0 /)
     aHCP(3,1:3) = (/ -0.5d0, -sqrt(3.d0)/2.d0, 0.d0 /)
     aHCP(4,1:3) = (/   0.d0,       0.d0, 1.5d0/coa  /)

     ! Normals of the 12 pyramidal planes (Miller indexes)
     nPy(1:4, 1) = (/  1,  0, -1,  1 /) ! Py1 fault
     nPy(1:4, 2) = (/  0,  1, -1,  1 /)
     nPy(1:4, 3) = (/ -1,  1,  0,  1 /)
     nPy(1:4, 4) = (/ -1,  0,  1,  1 /)
     nPy(1:4, 5) = (/  0, -1,  1,  1 /)
     nPy(1:4, 6) = (/  1, -1,  0,  1 /)
     nPy(1:4, 7) = (/  1,  0, -1, -1 /) ! Py1 fault
     nPy(1:4, 8) = (/  0,  1, -1, -1 /)
     nPy(1:4, 9) = (/ -1,  1,  0, -1 /)
     nPy(1:4,10) = (/ -1,  0,  1, -1 /)
     nPy(1:4,11) = (/  0, -1,  1, -1 /)
     nPy(1:4,12) = (/  1, -1,  0, -1 /)
     
     ! Patterns of the parent hcp lattice
     pattern0(1:3, 1,1) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
     pattern0(1:3, 2,1) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern0(1:3, 3,1) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern0(1:3, 4,1) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern0(1:3, 5,1) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern0(1:3, 6,1) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
     pattern0(1:3, 7,1) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern0(1:3, 8,1) = alat*(/   0.d0, -sqrt(3.d0)/3.d0,  0.5d0*coa /)
     pattern0(1:3, 9,1) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern0(1:3,10,1) = alat*(/  0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)
     pattern0(1:3,11,1) = alat*(/   0.d0, -sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern0(1:3,12,1) = alat*(/ -0.5d0,  sqrt(3.d0)/6.d0, -0.5d0*coa /)

     pattern0(1:3, 1,2) = alat*(/   1.d0,  0.d0           ,  0.d0 /)
     pattern0(1:3, 2,2) = alat*(/  0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern0(1:3, 3,2) = alat*(/  0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern0(1:3, 4,2) = alat*(/ -0.5d0,  sqrt(3.d0)/2.d0,  0.d0 /)
     pattern0(1:3, 5,2) = alat*(/ -0.5d0, -sqrt(3.d0)/2.d0,  0.d0 /)
     pattern0(1:3, 6,2) = alat*(/  -1.d0,  0.d0           ,  0.d0 /)
     pattern0(1:3, 7,2) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern0(1:3, 8,2) = alat*(/   0.d0,  sqrt(3.d0)/3.d0,  0.5d0*coa /)
     pattern0(1:3, 9,2) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0,  0.5d0*coa /)
     pattern0(1:3,10,2) = alat*(/  0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)
     pattern0(1:3,11,2) = alat*(/   0.d0,  sqrt(3.d0)/3.d0, -0.5d0*coa /)
     pattern0(1:3,12,2) = alat*(/ -0.5d0, -sqrt(3.d0)/6.d0, -0.5d0*coa /)

     ! Index of pyramidal planes
     IF (Present(nSys)) THEN
             IF ( (nSys.LT.1).OR.(nSys.GT.12) ) THEN
                     WRITE(0,'(a,i10)') 'nSys = ', nSys
                     STOP '< HcpPattern_AddPy1Twin >'
             END IF
             nSys1=nSys
             nSys2=nSys
     ELSE
             nSys1=1
             nSys2=12
     END IF

     ! Loop on pyramidal planes
     n1 = nPat -1
     n2 = nPat
     DO n=nSys1, nSys2
        ! Normal to pyramidal plane in cartesian coordinates
        normal(:) = MatMul( Dble(nPy(:,n)), aHCP(:,:) )
        inv_norm = 1.d0/Sqrt( Sum( normal(:)**2 ) )
        normal(:) = inv_norm*normal(:)
        ! New patterns indexes
        n1 = n1 + 2
        n2 = n2 + 2
        !   nickname
        labelPattern(n1) = "Twin hcp - plane A"
        labelPattern(n2) = "Twin hcp - plane B"
        !   number of neighbours
        nNeighPat(n1) = 12
        nNeighPat(n2) = 12
        !   loop on neighbours
        DO i=1, 12
           pattern(:,i,n1) = pattern0(:,i,1) - 2.d0*Sum( pattern0(:,i,1)*normal(:) )*normal(:)
           pattern(:,i,n2) = pattern0(:,i,2) - 2.d0*Sum( pattern0(:,i,2)*normal(:) )*normal(:)
        END DO
     END DO

     ! Reduced coordinates for each pattern
     DO n=nPat+1, nPat+2*(nSys2-nSys1+1)
        DO i=1, nNeighPat(n)
           uPattern(1:3,i,n) = pattern(1:3,i,n) / sqrt( Sum( pattern(1:3,i,n)**2) )
        END DO
     END DO

     ! Number of different atomic patterns in reference structure
     nPat = nPat+2*(nSys2-nSys1+1)

  END SUBROUTINE HcpPattern_AddPy1twin

  SUBROUTINE HcpPattern_AddPy1TB(alat,coa)
     ! Create patterns corresponding to 1st order pyramidal twin boundaries of hcp structure

     IMPLICIT NONE

     ! Lattice parameters and c/a ratio of hcp lattice
     REAL(kind(0.d0)), intent(in) :: alat, coa

     ! Rotation matrix
     REAL(kind(0.d0)), dimension(1:3,1:3) :: rot

     INTEGER :: i, n, n1, n2, n3, n4, nRot

     ! Allocation and initialization for pattern definition
     IF ( ( .NOT. Allocated(pattern) ).OR.( .NOT. Allocated(uPattern) ).OR.( nPat.LE.0 ) ) THEN
             WRITE(0,'(a)') 'You need to initialize HCP pattern before'
             STOP '< HcpPattern_AddPy1TB >'
     END IF

     IF (max_nNeighPat.LT.13) THEN
             WRITE(0,'(a,i0)') 'Problem with maximal number of neighbours in a pattern :  max_nNeighPat = ', max_nNeighPat
             WRITE(0,'(a)')    'It should be at least equal to 13'
             STOP '< HcpPattern_AddPy1TB >'
     END IF
     
     ! New patterns indexes
     n1 = nPat + 1
     n2 = nPat + 2
     n3 = nPat + 3
     n4 = nPat + 4

     ! Nickname
     labelPattern(n1) = "Twin boundary hcp - plane A - right"
     labelPattern(n2) = "Twin boundary hcp - plane B - right" 
     labelPattern(n3) = "Twin boundary hcp - plane A - left"
     labelPattern(n4) = "Twin boundary hcp - plane B - left" 

     ! Neighbour positions for pattern 1
     nNeighPat(n1) = 13
     pattern(:, 1,n1) = alat*(/ -0.5d0, -((Sqrt(3.d0)*(1.d0 + 2.d0*coa**2))/(6.d0 + 8.d0*coa**2)),     -((coa + 2.d0*coa**3)/(3.d0 + 4.d0*coa**2))  /)
     pattern(:, 2,n1) = alat*(/  0.5d0, -((Sqrt(3.d0)*(1.d0 + 2.d0*coa**2))/(6.d0 + 8.d0*coa**2)),     -((coa + 2.d0*coa**3)/(3.d0 + 4.d0*coa**2))  /)
     pattern(:, 3,n1) = alat*(/ -0.5d0, -(9.d0 + 13.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), coa/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 4,n1) = alat*(/  0.5d0, -(9.d0 + 13.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), coa/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 5,n1) = alat*(/ -0.5d0, -(3.d0 + 5.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)),  coa*(0.5d0 + 1.d0/(12.d0 + 16.d0*coa**2))  /)
     pattern(:, 6,n1) = alat*(/  0.5d0, -(3.d0 + 5.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)),  coa*(0.5d0 + 1.d0/(12.d0 + 16.d0*coa**2)) /)
     pattern(:, 7,n1) = alat*(/  0d0,   (Sqrt(3.d0)*(1.d0 + coa**2))/(3.d0 + 4.d0*coa**2),             (2.d0*(coa + coa**3))/(3.d0 + 4.d0*coa**2) /)
     pattern(:, 8,n1) = alat*(/ -0.5d0, (-3.d0 + 17.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (coa*(-15.d0 + 8.d0*coa**2))/(4.d0*(3.d0 + 4.d0*coa**2)) /)
     pattern(:, 9,n1) = alat*(/  0.5d0, (-3.d0 + 17.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (coa*(-15.d0 + 8.d0*coa**2))/(4.d0*(3.d0 + 4.d0*coa**2)) /)
     pattern(:,10,n1) = alat*(/ -0.5d0, (-9.d0 + 13.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (-25.d0*coa)/(12.d0 + 16.d0*coa**2)  /)
     pattern(:,11,n1) = alat*(/  0.5d0, (-9.d0 + 13.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (-25.d0*coa)/(12.d0 + 16.d0*coa**2)  /)
     pattern(:,12,n1) = alat*(/ -1.d0,  0.d0, 0.d0  /)
     pattern(:,13,n1) = alat*(/  1.d0,  0.d0, 0.d0  /)

     ! Neighbour positions for pattern 2
     nNeighPat(n2) = 11
     pattern(:, 1,n2) = alat*(/  0.d0,  -((Sqrt(3d0)*(1.d0 + coa**2))/(3.d0 + 4.d0*coa**2)),           (-2.d0*coa*(1.d0 + coa**2))/(3.d0 + 4.d0*coa**2)  /)
     pattern(:, 2,n2) = alat*(/ -0.5d0, -(9.d0 + 13.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), coa/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 3,n2) = alat*(/  0.5d0, -(9.d0 + 13.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), coa/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 4,n2) = alat*(/  0.d0,  -((Sqrt(3d0)*(2.d0 + 3.d0*coa**2))/(6.d0 + 8.d0*coa**2)),      coa*(0.5d0 + 1.d0/(12.d0 + 16.d0*coa**2))  /)
     pattern(:, 5,n2) = alat*(/ -0.5d0, (Sqrt(3.d0)*(1.d0 + 2.d0*coa**2))/(6.d0 + 8.d0*coa**2),        (coa + 2.d0*coa**3)/(3.d0 + 4.d0*coa**2)  /)
     pattern(:, 6,n2) = alat*(/  0.5d0, (Sqrt(3.d0)*(1.d0 + 2.d0*coa**2))/(6.d0 + 8.d0*coa**2),        (coa + 2.d0*coa**3)/(3.d0 + 4.d0*coa**2)  /)
     pattern(:, 7,n2) = alat*(/ -0.5d0, (-9.d0 + 13.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (-25.d0*coa)/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 8,n2) = alat*(/  0.5d0, (-9.d0 + 13.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (-25.d0*coa)/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 9,n2) = alat*(/  0.d0,  (Sqrt(3d0)*(-2.d0 + 7.d0*coa**2))/(6.d0 + 8.d0*coa**2),        (coa*(-23.d0 + 8.d0*coa**2))/(4.d0*(3.d0 + 4.d0*coa**2)) /)
     pattern(:,10,n2) = alat*(/ -1.d0, 0.d0, 0.d0  /)
     pattern(:,11,n2) = alat*(/  1.d0, 0.d0, 0.d0  /)

     ! Neighbour positioalat*ns for pattern 3
     nNeighPat(n3) = 11
     pattern(:, 1,n3) = alat*(/ -0.5d0, -((Sqrt(3.d0)*(1.d0 + 2.d0*coa**2))/(6.d0 + 8.d0*coa**2)),    -((coa + 2.d0*coa**3)/(3.d0 + 4.d0*coa**2))  /)
     pattern(:, 2,n3) = alat*(/  0.5d0, -((Sqrt(3.d0)*(1.d0 + 2.d0*coa**2))/(6.d0 + 8.d0*coa**2)),    -((coa + 2.d0*coa**3)/(3.d0 + 4.d0*coa**2)) /)
     pattern(:, 3,n3) = alat*(/  0.d0,   (6.d0 - 19.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (21.d0*coa - 8.d0*coa**3)/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 4,n3) = alat*(/ -0.5d0,  (9.d0 - 11.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (23.d0*coa)/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 5,n3) = alat*(/  0.5d0,  (9.d0 - 11.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), (23.d0*coa)/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 6,n3) = alat*(/  0.d0,   (Sqrt(3.d0)*(1.d0 + coa**2))/(3.d0 + 4.d0*coa**2),            (2.d0*(coa + coa**3))/(3.d0 + 4.d0*coa**2)  /)
     pattern(:, 7,n3) = alat*(/ -0.5d0,  (9.d0 + 11.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), coa/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 8,n3) = alat*(/  0.5d0,  (9.d0 + 11.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), coa/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 9,n3) = alat*(/  0.d0,   (6.d0 +  7.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)), coa*(-0.5d0 + 1.d0/(12.d0 + 16.d0*coa**2)) /)
     pattern(:,10,n3) = alat*(/ -1.d0, 0.d0, 0.d0  /)
     pattern(:,11,n3) = alat*(/  1.d0, 0.d0, 0.d0  /)

     ! Neighbour positioalat*ns for pattern 4
     nNeighPat(n4) = 13
     pattern(:, 1,n4) = alat*(/  0.d0, -Sqrt(3.d0)*(1.d0 +      coa**2)/(3.d0 + 4.d0*coa**2),(-2.d0*coa*(1.d0 + coa**2))/(3.d0 + 4.d0*coa**2) /)
     pattern(:, 2,n4) = alat*(/ -0.5d0, Sqrt(3.d0)*(1.d0 - 5.d0*coa**2)/(6.d0 + 8.d0*coa**2),(13.d0*coa - 8.d0*coa**3)/(12.d0 + 16.d0*coa**2) /)
     pattern(:, 3,n4) = alat*(/  0.5d0, Sqrt(3.d0)*(1.d0 - 5.d0*coa**2)/(6.d0 + 8.d0*coa**2),(13.d0*coa - 8.d0*coa**3)/(12.d0 + 16.d0*coa**2) /)
     pattern(:, 4,n4) = alat*(/ -0.5d0, (9.d0 - 11.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)),(23.d0*coa)/(12.d0 + 16.d0*coa**2) /)
     pattern(:, 5,n4) = alat*(/  0.5d0, (9.d0 - 11.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)),(23.d0*coa)/(12.d0 + 16.d0*coa**2) /)
     pattern(:, 6,n4) = alat*(/ -0.5d0, (Sqrt(3.d0)*(1.d0 + 2.d0*coa**2))/(6.d0 + 8.d0*coa**2),(coa + 2.d0*coa**3)/(3.d0 + 4.d0*coa**2) /)
     pattern(:, 7,n4) = alat*(/  0.5d0, (Sqrt(3.d0)*(1.d0 + 2.d0*coa**2))/(6.d0 + 8.d0*coa**2),(coa + 2.d0*coa**3)/(3.d0 + 4.d0*coa**2) /)
     pattern(:, 8,n4) = alat*(/ -0.5d0, (9.d0 + 11.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)),coa/(12.d0 + 16.d0*coa**2)  /)
     pattern(:, 9,n4) = alat*(/  0.5d0, (9.d0 + 11.d0*coa**2)/(2.d0*Sqrt(3.d0)*(3.d0 + 4.d0*coa**2)),coa/(12.d0 + 16.d0*coa**2)  /)
     pattern(:,10,n4) = alat*(/ -0.5d0, (Sqrt(3.d0)*(1.d0 + coa**2))/(6.d0 + 8.d0*coa**2),coa*(-0.5d0 + 1/(12.d0 + 16.d0*coa**2)) /)
     pattern(:,11,n4) = alat*(/  0.5d0, (Sqrt(3.d0)*(1.d0 + coa**2))/(6.d0 + 8.d0*coa**2),coa*(-0.5d0 + 1/(12.d0 + 16.d0*coa**2)) /)
     pattern(:,12,n4) = alat*(/ -1.d0, 0.d0, 0.d0   /)
     pattern(:,13,n4) = alat*(/  1.d0, 0.d0, 0.d0   /)

     ! Rotation pi/3
     rot(1:3,1:3) = Reshape( &
             (/ 0.5d0, -0.866025403784439d0, 0.d0, &
                0.866025403784439d0,  0.5d0, 0.d0, &
                0.d0, 0.d0, 1.d0 /), (/3,3/) )
     n = n4
     DO nRot=1, 5*4
        n = n + 1
        labelPattern(n) = labelPattern(n-4)
        nNeighPat(n) = nNeighPat(n-4)
        DO i=1, nNeighPat(n)
           pattern(:,i,n) = MatMul( rot(:,:), pattern(:,i,n-4) )
        END DO
     END DO

     !! ! Mirror in basal plane
     !! rot(1:3,1:3) = Reshape( &
     !!         (/ 1.d0, 0.d0, 0.d0, &
     !!            0.d0, 1.d0, 0.d0, &
     !!            0.d0, 0.d0,-1.d0 /), (/3,3/) )
     !! DO n=n1, n1+23
     !!    labelPattern(n+24) = labelPattern(n)
     !!    nNeighPat(n+24) = nNeighPat(n)
     !!    DO i=1, nNeighPat(n)
     !!       pattern(:,i,n+24) = MatMul( rot(:,:), pattern(:,i,n) )
     !!    END DO
     !! END DO
        

     ! Reduced coordinates for each pattern
     DO n=nPat+1, nPat+24
        DO i=1, nNeighPat(n)
           uPattern(1:3,i,n) = pattern(1:3,i,n) / sqrt( Sum( pattern(1:3,i,n)**2) )
        END DO
     END DO

     ! Number of different atomic patterns in reference structure
     nPat = nPat+24

  END SUBROUTINE HcpPattern_AddPy1TB

  SUBROUTINE InitGradElasticDisplacement(imm)

     IMPLICIT NONE

     ! Maximal number of atoms
     INTEGER, intent(in) :: imm

     IF (Allocated(gradElasticDisplacement)) Deallocate(gradElasticDisplacement)
     ALLOCATE(gradElasticDisplacement(1:3,1:3,1:imm))
     gradElasticDisplacement(:,:,:) = 0.d0

  END SUBROUTINE InitGradElasticDisplacement

  SUBROUTINE BuildGradElasticDisplacement(im, xp, at, nNeigh, iNeigh, out)
     ! Calculate elastic displacement gradient from finite differences 
     ! by comparing input structure with registered patterns
     
     USE Babel_data, ONLY : verbosity
     IMPLICIT NONE

     ! Number of atoms
     INTEGER, intent(in) :: im
     ! Atoms coordinates in input configuration
     REAL(kind(0.d0)), dimension(:,:), intent(in) :: xp
     ! Periodicity vectors of the input structure: at(1:3,1), at(1:3,2), ...
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
     ! Number of neighbours for each atom 
     INTEGER, dimension(:), intent(in) :: nNeigh
     ! Neighbour indexes for each atom i: iNeigh(1:nNeigh(i), i)
     INTEGER, dimension(:,:), intent(in) :: iNeigh
     ! Output unit 
     INTEGER, intent(in), optional :: out

     ! Local variables
     INTEGER :: nNeighMax, nMax, npMax, i, j, n, n0, p, p0, np, np0, nn, nUnmatched, nUnmatched0, err
     REAL(kind(0.d0)) :: f, f0, dFtotal2min, dFtotal2, dR2, dR2min
     INTEGER, dimension(1:2) :: iPos
     INTEGER, dimension(1:max_nPat) :: nMatch
     REAL(kind(0.d0)) :: invR
     REAL(kind(0.d0)), dimension(1:3) :: ds
     REAL(kind(0.d0)), dimension(:,:), allocatable :: R, u
     REAL(kind(0.d0)), dimension(:,:), allocatable :: cosTab
     REAL(kind(0.d0)), dimension(:,:,:), allocatable :: Rordered, R0ordered
     REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
     REAL(kind(0.d0)), dimension(1:3,1:3) :: Amat, Bmat, invBmat, Ftotal
     ! Threshold cosine for pairing vectors
     REAL(kind(0.d0)) :: cosThreshold 
     ! Number of atoms in pattern p: nAtoms_in_pat(p)
     INTEGER, dimension(1:max_nPat) :: nAtoms_in_pat


     ! Threshold cosine for pairing vectors
     cosThreshold = cos(patternAngleThreshold/180.d0*3.14159265358979d0)

     ! Check number of registered patterns
     IF (nPat.LE.0) THEN
             WRITE(0,'(a)') 'Patterns need to be identified first in reference configuration'
             STOP '< BuildGradDisplacement >'
     END IF

     ! Inverse periodicity vectors
     IF (perioElasticDisplacement) THEN
             CALL Mat3Inv(at, inv_at)
     ELSE
             inv_at(:,:)=0.d0
     END IF

     ! Maximal number of neighbours in input structure
     nNeighMax = MaxVal(nNeigh(1:im))
     ! Matrix containing neighbour vectors (and normalized vectors)
     ALLOCATE(R(1:3, 1:nNeighMax))
     ALLOCATE(u(1:3, 1:nNeighMax))

     ! Maximal number of neighbours in patterns
     npMax = MaxVal(nNeighPat(1:nPat))

     ! Table containing cosine of neighbour vectors between input and reference
     ! structures
     ALLOCATE( cosTab(1:nNeighMax, 1:npMax) )

     ! Matrix containing ordered list of vectors in input and reference
     ! structures for each pattern
     ALLOCATE( Rordered(1:3, 1:nNeighMax, 1:nPat) )
     ALLOCATE( R0ordered(1:3, 1:nNeighMax, 1:nPat) )

     !-- Initialization -----------------
     ! table recording pattern number for atom n in reference structure 
     IF (Allocated(inpPattern)) Deallocate(inpPattern)
     Allocate(inpPattern(1:im))
     inpPattern(:)=0
     ! table counting number of atoms in each pattern
     nAtoms_in_pat(:) = 0

     !-------------------------------------------
     ! Loop on atoms
     atom_loop: DO i=1, im

        ! Number of neighbours for current atom
        nMax = nNeigh(i)

        IF (nMax.LE.0) THEN
                WRITE(0,'(a,i0,a)') 'Atom ', i, ' does not have neighbours'
                STOP '< BuildGradDisplacement >'
        END IF

        !----------------------------------------
        ! Build neighbour tables for current atom
        neighbour_loop: DO n=1, nMax
          ! Index of neighbour atom
          j = iNeigh(n, i)
          R(1:3,n) = xp(1:3,j) - xp(1:3,i)
          IF (perioElasticDisplacement) THEN
                  ds(:) = MatMul( inv_at(:,:), R(1:3,n) )
                  ds(:) = aNInt(ds(:))
                  R(1:3,n) = R(1:3,n) - MatMul( at(:,:), ds(:) )
          END IF
          ! Normalize neighbour vector
          invR = 1.d0/sqrt( Sum( R(1:3,n)**2 ) )
          u(1:3,n) = R(1:3,n)*invR
        END DO neighbour_loop

        !----------------------------------------
        ! Look for closest pattern
        pattern_loop: DO p=1, nPat

           ! Number of neigbours for current pattern
           npMax = nNeighPat(p)

           ! Table containing cosine of neighbour vectors 
           ! between input and reference structures
           cosTab(:,:) = -1.d0
           DO n=1, nMax
              DO np=1, npMax
                 cosTab(n,np) = Sum( u(:,n)*uPattern(:,np,p) )
              END DO
           END DO

           ! Look for match between vectors in input and reference config.
           n = 0
           DO 
              ! Index where the maximal cosinus is obtained
              iPos(:) = MaxLoc( cosTab(1:nMax, 1:npMax) )
              n0 = iPos(1) 
              np0 = iPos(2)

              ! Check if the two vectors are not too different
              IF (cosTab(n0,nP0).LT.cosThreshold) EXIT
              
              ! We found a new match
              n = n + 1
              Rordered(1:3,n,p) = R(1:3,n0)
              R0ordered(1:3,n,p) = pattern(1:3,np0,p)

              ! Check if the number of vectors in input structure has not been exhausted
              IF (n.GE.nMax) EXIT
              ! ... same for reference structure
              IF (n.GE.npMax) EXIT

              ! Cancel corresponding line and column in cosine table
              cosTab(n0,:) = -1.d0
              cosTab(:,np0) = -1.d0

           END DO
           ! Number of vectors that match between input structure and pattern p
           nMatch(p) = n

        END DO pattern_loop

        SELECT CASE(patternSelectionMethod)
        CASE(1)
                ! Select pattern which maximises the number of matchs 
                !   look to the NORM OF THE DISPLACEMENT if 2 patterns are equivalent
                n0 = 0    !  Number of match for best pattern
                f0 = 0.d0 ! Best pattern macth value
                p0 = 0    ! Best pattern index
                DO p=1, nPat
                   n = nMatch(p)
                   f = dble(n)/dble(nNeighPat(p))
                   dR2 = Sum( ( Rordered(1:3,1:n,p) - R0ordered(1:3,1:n,p) )**2/dble(n) )
                   IF (f.GT.f0) THEN
                           n0 = n
                           f0 = f
                           p0 = p
                           dR2min = dR2
                   ELSE IF (f.EQ.f0) THEN
                           IF ( dR2 .LT. dR2min ) THEN
                                   n0 = n
                                   f0 = f
                                   p0 = p
                                   dR2min = dR2
                           END IF
                   END IF
                   IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) &                  ! DEBUG
                           WRITE(out,'(4(a,i3),a,g13.6)') 'atom ', i, ', pattern ', p, ' - ', & ! DEBUG
                                  n, '/', nNeighPat(p), ' matchs - dR2 = ', dR2                 ! DEBUG
                END DO

                ! Atom i belongs to pattern p0
                inpPattern(i) = p0
                IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) THEN  ! DEBUG
                        WRITE(out,'(a,i0)') '  Selected pattern: ', p0                  ! DEBUG
                        WRITE(out,*)                                                    ! DEBUG
                END IF
        CASE(2)
                ! Select pattern which maximises the number of matchs 
                !   look to the norm of the ELASTIC GRADIENT if 2 patterns are equivalent
                n0 = 0    ! Number of match for best pattern
                f0 = 0.d0 ! Best pattern macth value
                p0 = 0    ! Best pattern index
                DO p=1, nPat
                   n=nMatch(p)
                   f = dble(n)/dble(nNeighPat(p))
                   Amat(:,:) = MatMul( R0ordered(1:3,1:n,p), Transpose( Rordered(1:3,1:n,p) - R0ordered(1:3,1:n,p) ) )
                   Bmat(:,:) = MatMul( R0ordered(1:3,1:n,p), Transpose( R0ordered(1:3,1:n,p) ) )
                   CALL Mat3Inv(Bmat, invBmat, err)
                   IF (err.NE.0) CYCLE
                   Ftotal(:,:) = MatMul( Amat(:,:), invBmat(:,:) )
                   dFtotal2 = Sum( Ftotal(:,:)**2 )
                   IF (f.GT.f0) THEN
                           n0 = n
                           f0 = f
                           p0 = p
                           dFtotal2min = dFtotal2
                   ELSE IF (f.EQ.f0) THEN
                           IF ( dFtotal2.LT.dFtotal2min ) THEN
                                   n0 = n
                                   f0 = f
                                   p0 = p
                                   dFtotal2min = dFtotal2
                           END IF
                   END IF
                   IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) &                  ! DEBUG
                          WRITE(out,'(4(a,i3),a,g13.6)') 'atom ', i, ', pattern ', p, ' - ', &
                                   n, '/', nNeighPat(p), ' matchs - |dF|^2 = ', dFtotal2 !DEBUG
                END DO

                ! Atom i belongs to pattern p0
                inpPattern(i) = p0
                IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) THEN  ! DEBUG
                        WRITE(out,'(a,i0)') '  Selected pattern: ', p0                  ! DEBUG
                        WRITE(out,*)                                                    ! DEBUG
                END IF
        CASE(3)
                ! Select pattern which minimizes the number of unmatched bonds 
                !   look to the NORM OF THE DISPLACEMENT if 2 patterns are equivalent
                n0 = 0    !  Number of match for best pattern
                nUnmatched0 = MaxVal( nNeighPat(1:nPat) ) ! Number of unmatched bonds
                p0 = 0    ! Best pattern index
                DO p=1, nPat
                   n = nMatch(p)
                   nUnmatched = nNeighPat(p) - n
                   dR2 = Sum( ( Rordered(1:3,1:n,p) - R0ordered(1:3,1:n,p) )**2/dble(n) )
                   IF (nUnmatched.LT.nUnmatched0) THEN
                           n0 = n
                           nUnmatched0 = nUnmatched
                           p0 = p
                           dR2min = dR2
                   ELSE IF (nUnmatched.EQ.nUnmatched0) THEN
                           IF ( dR2 .LT. dR2min ) THEN
                                   n0 = n
                                   nUnmatched0 = nUnmatched
                                   p0 = p
                                   dR2min = dR2
                           END IF
                   END IF
                   IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) &                  ! DEBUG
                           WRITE(out,'(4(a,i3),a,g13.6)') 'atom ', i, ', pattern ', p, ' - ', & ! DEBUG
                                  n, '/', nNeighPat(p), ' matchs - dR2 = ', dR2                 ! DEBUG
                END DO

                ! Atom i belongs to pattern p0
                inpPattern(i) = p0
                IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) THEN          ! DEBUG
                        WRITE(out,'(a,i0)') '  Selected pattern: ', p0                  ! DEBUG
                        WRITE(out,*)                                                    ! DEBUG
                END IF
        CASE(4)
                ! Select pattern which minimizes the number of unmatched bonds 
                !   look to the norm of the ELASTIC GRADIENT if 2 patterns are equivalent
                n0 = 0    !  Number of match for best pattern
                nUnmatched0 = MaxVal( nNeighPat(1:nPat) ) ! Number of unmatched bonds
                p0 = 0    ! Best pattern index
                DO p=1, nPat
                   n = nMatch(p)
                   nUnmatched = nNeighPat(p) - n
                   IF (nUnmatched.GT.nUnmatched0) CYCLE
                   Amat(:,:) = MatMul( R0ordered(1:3,1:n,p), Transpose( Rordered(1:3,1:n,p) - R0ordered(1:3,1:n,p) ) )
                   Bmat(:,:) = MatMul( R0ordered(1:3,1:n,p), Transpose( R0ordered(1:3,1:n,p) ) )
                   CALL Mat3Inv(Bmat, invBmat, err)
                   IF (err.NE.0) CYCLE
                   Ftotal(:,:) = MatMul( Amat(:,:), invBmat(:,:) )
                   dFtotal2 = Sum( Ftotal(:,:)**2 )
                   IF (nUnmatched.LT.nUnmatched0) THEN
                           n0 = n
                           nUnmatched0 = nUnmatched
                           p0 = p
                           dFtotal2min = dFtotal2
                   ELSE IF (nUnmatched.EQ.nUnmatched0) THEN
                           IF ( dFtotal2.LT.dFtotal2min ) THEN
                                   n0 = n
                                   nUnmatched0 = nUnmatched
                                   p0 = p
                                   dFtotal2min = dFtotal2
                           END IF
                   END IF
                   IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) &                  ! DEBUG
                          WRITE(out,'(4(a,i3),a,g13.6)') 'atom ', i, ', pattern ', p, ' - ', &
                                   n, '/', nNeighPat(p), ' matchs - |dF|^2 = ', dFtotal2 !DEBUG
                END DO

                ! Atom i belongs to pattern p0
                inpPattern(i) = p0
                IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) THEN          ! DEBUG
                        WRITE(out,'(a,i0)') '  Selected pattern: ', p0                  ! DEBUG
                        WRITE(out,*)                                                    ! DEBUG
                END IF
        END SELECT

        ! Pattern p0 counts one more atom
        nAtoms_in_pat(p0) = nAtoms_in_pat(p0) + 1

        !--------------------------------------------------------------------
        ! Calculate gradient of displacement
        Amat(:,:) = MatMul( R0ordered(1:3,1:n0,p0), Transpose( Rordered( 1:3,1:n0,p0) - R0ordered(1:3,1:n0,p0) ) )
        Bmat(:,:) = MatMul( R0ordered(1:3,1:n0,p0), Transpose( R0ordered(1:3,1:n0,p0) ) )
        CALL Mat3Inv(Bmat, invBmat)
        Ftotal(:,:) = MatMul( Amat(:,:), invBmat(:,:) )
        gradElasticDisplacement(1:3,1:3,i) = Transpose( Ftotal(:,:) )

        IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) THEN
                WRITE(out,'(4(a,i0),a)') ' atom ', i, ' corresponds to pattern ', &
                        p0, ' (', n0, ' matchs / ', nNeigh(i), ' neighbours)'
                DO nn=1, n0
                   WRITE(out,'(a,i3,2(a,3g14.6))') '  neighbour ', nn, &
                           ' R(1:3) = ',  Rordered(1:3,nn,p0), &
                           ' - R0(1:3) = ', R0ordered(1:3,nn,p0)
                END DO
                WRITE(out,'(a,3g14.6)') '  dUx / dR(1:3) = ', Ftotal(1,1:3)
                WRITE(out,'(a,3g14.6)') '  dUy / dR(1:3) = ', Ftotal(2,1:3)
                WRITE(out,'(a,3g14.6)') '  dUz / dR(1:3) = ', Ftotal(3,1:3)
                WRITE(out,*)
        END IF

     END DO atom_loop

     DEALLOCATE(R) ; DEALLOCATE(u) ; DEALLOCATE(cosTab)
     DEALLOCATE(Rordered) ; DEALLOCATE(R0ordered)

     IF (verbosity.GE.verbosity_max) THEN
             WRITE(out,*)
             DO p=1, nPat
                WRITE(out,'(2(a,i0))') 'Number of atoms corresponding to pattern ', p,': ', &
                        nAtoms_in_pat(p)
             END DO
             WRITE(out,*)
     END IF

  END SUBROUTINE BuildGradElasticDisplacement

  SUBROUTINE Mat3Inv(A, B, err)
    ! return matrix B which is the inverse of matrix A
    ! A and B are 3x3 matrices
    ! err=0  if inversion is possible
    !     -1 if the determinant is equal to 0

    implicit none

    REAL(kind(0.d0)), dimension(3,3), intent(in) :: A
    REAL(kind(0.d0)), dimension(3,3), intent(out) :: B
    INTEGER, intent(out), optional :: err

    REAL(kind(0.d0)) :: invdet, det, tol
    REAL(kind(0.d0)), parameter :: zero=1.d-12

    ! Tolerance for zero determinant
    tol = zero*MaxVal( Abs( a(:,:) ) )
       
    b(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
    b(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
    b(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
       
    b(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
    b(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
    b(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)
       
    b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
    b(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    det =  a(1,1)*b(1,1) + a(1,2)*b(2,1) + a(1,3)*b(3,1) 
    IF (Abs(det).GT.tol) THEN
            invdet = 1.d0/det
            b(:,:)=b(:,:)*invdet
            IF (Present(err)) err=0
    ELSE
            b(:,:)=0.d0
            IF (Present(err)) err=-1
    END IF

  END SUBROUTINE mat3inv

END MODULE GradElasticDisplacementModule

