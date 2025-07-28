MODULE GradElasticDisplacementModule

   IMPLICIT NONE

   SAVE
   
   ! If patternRead is set to true, the pattern is read in the file patternFile
   LOGICAL :: patternRead
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


   ! Number of different atomic patterns in reference structure (and maximal ! value)
   INTEGER, parameter, private :: max_nPat=20
   INTEGER, private :: nPat=0

   ! Atomic patterns definition:
   !   nNeighPat(n): number of neighbours in pattern n
   INTEGER, dimension(1:max_nPat), private :: nNeighPat
   !   pattern(1:3,i,n): coordinate of neighbour i for pattern n
   REAL(kind(0.d0)), dimension(:,:,:), allocatable, private :: pattern
   !   corresponding normalized vectors
   REAL(kind(0.d0)), dimension(:,:,:), allocatable, private :: uPattern

   ! refPattern(n) and inpPattern(n): pattern number for atom n in reference and input structures
   INTEGER, dimension(:), allocatable, private :: refPattern
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

     USE Babel_data, ONLY : verbosity
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
     INTEGER :: nNeighMax, i, n, j, nMax, np
     LOGICAL :: new_pattern
     REAL(kind(0.d0)) :: invR0
     REAL(kind(0.d0)), dimension(1:3) :: ds0
     REAL(kind(0.d0)), dimension(:,:), allocatable :: R0
     REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at0

     IF (verbosity.GE.verbosity_max) THEN
             WRITE(out,*)
             WRITE(out,'(a)') 'Look for pattern definition in reference structure'
             WRITE(out,*)
     END IF

     !-----------------------------------
     ! table recording pattern number for atom n in reference structure 
     IF (Allocated(refPattern)) Deallocate(refPattern)
     Allocate(refPattern(1:im0))
     refPattern(:)=0

     ! Maximal number of neighbours
     nNeighMax = MaxVal(nNeigh(1:im0))

     ! Allocation and initialization for pattern definition
     nPat=0
     IF (Allocated(pattern)) Deallocate(pattern)
     ALLOCATE(pattern(1:3,1:nNeighMax,1:max_nPat))
     IF (Allocated(uPattern)) Deallocate(uPattern)
     ALLOCATE(uPattern(1:3,1:nNeighMax,1:max_nPat))
     pattern(:,:,:)=0.d0
     uPattern(:,:,:)=0.d0
     nNeighPat(:)=0

     ! Local array
     ALLOCATE(R0(1:3,1:nNeighMax))

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

                IF (verbosity.GE.verbosity_debug_pattern) THEN
                        WRITE(out,*)
                        WRITE(out,'(a)') 'New pattern found'
                        CALL PrintPattern(nPat,out)
                        WRITE(out,*)
                END IF
        ELSE 
                ! Atom i corresponds to pattern np
                refPattern(i) = np

                IF (verbosity.GE.verbosity_debug_pattern) &
                        WRITE(out,'(a,i0)') 'Pattern already exists: equivalent to pattern ', np
        END IF

     END DO atom_loop

     Deallocate(R0)

     IF (verbosity.GE.verbosity_max) THEN
             CALL PrintAllPatterns(out)
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
    WRITE(out,'(a,i0)') '  number of neighbours in pattern definition: ', nNeighPat(n)
    DO i=1, nNeighPat(n)
       WRITE(out,'(a,i3,a,3g14.6)')'    atom ', i, ': R(1:3) = ', pattern(1:3,i,n)
    END DO

  END SUBROUTINE PrintPattern

  SUBROUTINE ReadPattern(inp, out)
     ! Read pattern definition in file connected to unit inp

     USE Babel_data, ONLY : verbosity
     IMPLICIT NONE
     ! input and output units
     INTEGER, intent(in) :: inp, out

     INTEGER :: i, n, nNeigh, nNeighMax

     ! Number of different atomic patterns in reference structure
     CALL Comment(inp)
     READ(inp,*) nPat
     IF ( (nPat.LE.0) .OR. (nPat.GT.max_nPat) ) THEN
             WRITE(0,'(a,i0)') 'nPat = ', nPat
             STOP '< ReadPattern >'
     END IF

     ! First pass to determine the size needed for the arrays
     nNeighMax=0
     DO n=1, nPat
        ! number of neighbours in pattern n
        CALL Comment(inp)
        READ(inp,*) nNeigh
        IF (nNeigh.GT.nNeighMax) nNeighMax=nNeigh
        ! coordinate of neighbour i for pattern n
        DO i=1, nNeigh
           CALL Comment(inp)
           READ(inp,*) 
        END DO
     END DO
        
     REWIND(inp)

     ! Allocation and initialization for pattern definition
     IF (Allocated(pattern)) Deallocate(pattern)
     ALLOCATE(pattern(1:3,1:nNeighMax,1:max_nPat))
     IF (Allocated(uPattern)) Deallocate(uPattern)
     ALLOCATE(uPattern(1:3,1:nNeighMax,1:max_nPat))
     pattern(:,:,:)=0.d0
     uPattern(:,:,:)=0.d0
     nNeighPat(:)=0

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
     INTEGER :: nNeighMax, nMax, npMax, i, j, n, n0, p, p0, np, np0, nn
     INTEGER, dimension(1:2) :: iPos
     INTEGER, dimension(1:max_nPat) :: nMatch
     REAL(kind(0.d0)) :: invR, dR2min, dR2
     REAL(kind(0.d0)), dimension(1:3) :: ds
     REAL(kind(0.d0)), dimension(:,:), allocatable :: R, u
     REAL(kind(0.d0)), dimension(:,:), allocatable :: cosTab
     REAL(kind(0.d0)), dimension(:,:,:), allocatable :: Rordered, R0ordered
     REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
     REAL(kind(0.d0)), dimension(1:3,1:3) :: Amat, Bmat, invBmat, Ftotal
     ! Threshold cosine for pairing vectors
     REAL(kind(0.d0)) :: cosThreshold 


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

     !-----------------------------------
     ! table recording pattern number for atom n in reference structure 
     IF (Allocated(inpPattern)) Deallocate(inpPattern)
     Allocate(inpPattern(1:im))
     inpPattern(:)=0
     !-------------------------------------------
     ! Loop on atoms
     atom_loop: DO i=1, im

        ! Number of neighbours for current atom
        nMax = nNeigh(i)

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

        !!$ ! Select pattern which minimises the strain
        !!$ n0 = 0   ! Number of match for best pattern
        !!$ p0 = 0  ! Best pattern index
        !!$ dR2min = 0.d0
        !!$ DO p=1, nPat
        !!$    n = nMatch(p)
        !!$    dR2 = Sum( ( Rordered(1:3,1:n,p) - R0ordered(1:3,1:n,p))**2)/dble(n**2)
        !!$    IF (p.EQ.1) THEN
        !!$            n0 = n
        !!$            p0 = p
        !!$            dR2min = dR2
        !!$    ELSE IF (dR2.LT.dR2min) THEN
        !!$            n0 = n
        !!$            p0 = p
        !!$            dR2min = dR2
        !!$    END IF
        !!$ END DO

        ! Select pattern which maximises the number of matchs (tolerance +/-1)
        n0 = 0   ! Number of match for best pattern
        p0 = 0  ! Best pattern index
        DO p=1, nPat
           n=nMatch(p)
           IF (n.GT.n0+0) THEN
                   n0 = n
                   p0 = p
           ELSE IF (n.GE.n0-0) THEN
                   IF ( Sum( ( Rordered(1:3,1:n,p) - R0ordered(1:3,1:n,p) )**2 ) &
                           .LT.Sum( ( Rordered(1:3,1:n0,p0) - R0ordered(1:3,1:n0,p0) )**2 ) ) &
                           n0 = n
                           p0 = p
           END IF
        END DO

        inpPattern(i) = p0
        IF (verbosity.GE.verbosity_max) &
                WRITE(out,'(4(a,i0),a)') ' atom ', i, ' corresponds to pattern ', &
                        p0, ' (', n0, ' matchs / ', nNeigh(i), ' neighbours)'

        ! Calculate gradient of displacement
        Amat(:,:) = MatMul( R0ordered(1:3,1:n0,p0), Transpose( Rordered(1:3,1:n0,p0) - R0ordered(1:3,1:n0,p0) ) )
        Bmat(:,:) = MatMul( R0ordered(1:3,1:n0,p0), Transpose( R0ordered(1:3,1:n0,p0) ) )
        CALL Mat3Inv(Bmat, invBmat)
        Ftotal(:,:) = MatMul( Amat(:,:), invBmat(:,:) )
        gradElasticDisplacement(1:3,1:3,i) = Transpose( Ftotal(:,:) )

        IF (verbosity.GE.verbosity_debug_gradElasticDisplacement) THEN
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

  END SUBROUTINE BuildGradElasticDisplacement

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

END MODULE GradElasticDisplacementModule

