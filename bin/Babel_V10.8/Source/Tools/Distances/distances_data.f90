MODULE Distances_module

  SAVE

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter :: Distance_Zero2=Distance_Zero*Distance_Zero

  !========================================================
  ! Radius of the neighbourhood sphere
  REAL(kind(0.d0)) :: rNeigh

  ! Index of the atom for which we calculate the neighbourhood
  INTEGER :: iAtom

  !========================================================
  ! If duplicate=.true. the unit cell is duplicated lat(1), lat(2)... times in each
  ! direction given by at(1:3,1), at(1:3,2)...
  LOGICAL :: duplicate
  INTEGER, dimension(3) :: lat
  
  ! Name of the input file defining the dislocation
  CHARACTER(len=100) :: input_file

  !========================================================
  ! DEFINITION OF THE STRUCTURE
  ! Lattice vector coordinates (A): at(1:3,1), at(2:3,2), ...
  ! (input and reference structures)
  LOGICAL :: at_defined
  REAL(kind(0.d0)), dimension(1:3, 1:3) :: at
  ! Maximal number of atoms in simulation box (used to allocate tables)
  INTEGER :: imm=0
  ! Real number of atoms in simulation box
  INTEGER :: im=0
  ! Atom real coordinates: xp(1:3,:) 
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp
  ! Atom type 
  INTEGER, dimension(:), allocatable :: iTyp

  INTEGER :: verbosity
  INTEGER, parameter :: verbosity_debug=10
  LOGICAL :: debug

CONTAINS

  SUBROUTINE Read_distances(out)

     USE structure_module
     IMPLICIT NONE
     INTEGER, intent(in) :: out

     ! Format of the input structure files
     LOGICAL :: inpXyz, inpCfg, inpGin, inpLisa, inpSiesta, inpNDM
     ! Input structure file
     CHARACTER(len=100) :: inpFile
     ! Type of input file
     CHARACTER(len=6) :: fileType

     REAL(kind(0.d0)) :: alat, at_norm2
     REAL(kind(0.d0)), dimension(1:3,1:3) :: at_temp


     NAMELIST /input/ alat, at, &
             imm, duplicate, lat, &
             inpXyz, inpCfg, inpGin, inpLisa, inpSiesta, inpNDM, inpFile, &
             rNeigh, iAtom, &
             verbosity, debug

     alat=1.d0
     at(1:3,1:3)=0.d0     
     imm=0
     duplicate=.FALSE.
     inpXyz=.FALSE. ; inpCfg=.FALSE. ; inpGin=.FALSE. ; inpLisa=.FALSE. ; inpSiesta=.FALSE. ; inpNDM=.FALSE.
     lat(1:3)=1
     rNeigh=-1.d0
     iAtom=0
     verbosity=4
     debug=.FALSE.

    !===========================================
    ! Read input data
    OPEN(file=input_file,unit=50,action='read',status='old')
    READ(50,nml=input)
    CLOSE(50)

    !===========================================
    ! Basis vectors could be defined here or read in structure file
    at_norm2=SUM( at(1:3,1:3)**2 )
    IF (at_norm2.GT.Distance_Zero2) THEN
            ! Basis vectors defined
            at(:,:) = alat*at(:,:)
            at_defined=.TRUE.
    ELSE
            at_defined=.FALSE.
    END IF

    !===========================================
    ! Read input structure file

    IF ( Count( (/inpXyz,inpCfg,inpGin,inpLisa,inpSiesta,inpNDM/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin, Lisa, and Siesta formats for input structure'
            STOP '< Read_Distance >'
    END IF

    IF (inpXyz.OR.inpCfg.OR.inpGin.OR.inpLisa.OR.inpSiesta.OR.inpNDM) THEN
            IF (inpXyz) fileType='xyz'
            IF (inpCfg) fileType='cfg'
            IF (inpGin) fileType='gin'
            IF (inpSiesta) fileType='siesta'
            IF (inpLisa) fileType='lisa'
            IF (inpNDM) fileType='ndm'
            CALL ReadStructure(xp, iTyp, im, imm, at_temp, inpFile, fileType)
    END IF

    !===========================================
    ! Basis vectors could be defined here or read in structure file
    at_norm2=SUM( at_temp(1:3,1:3)**2 )
    IF (at_norm2.GT.Distance_Zero2) THEN
            IF (at_defined) THEN
                    WRITE(0,'(a)') 'Periodicity vectors at(1:3,i)&
                        & defined in input and structure files '
                    STOP '< Read_Displacement >'
            ELSE
                    at(:,:) = at_temp
                    at_defined=.TRUE.
            END IF
    END IF


  END SUBROUTINE Read_distances

  SUBROUTINE Calculate_Distances(dNeigh, nNeigh)
      ! Store in table dNeigh(:) the distance between atom iAtom
      ! and all other atom in the simulation box if the distance 
      ! is smaller than rNeigh

      USE structure_module
      
      IMPLICIT NONE
      
      REAL(kind(0.d0)), dimension(:), intent(out) :: dNeigh
      INTEGER, intent(out) :: nNeigh

      INTEGER, dimension(1) :: nMinVec
      INTEGER :: i, n, nMin
      REAL(kind(0.d0)), dimension(1:3) :: xp0, dxp, dxc
      REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
      REAL(kind(0.d0)) :: at_norm2, rNeigh2, d2, d0, dMin

      IF ( (iAtom.LE.0).OR.(iAtom.GT.im) ) THEN
              WRITE(0,'(a,i0)') 'iAtom = ', iAtom
              STOP "< Calculate_Distances >"
      END IF

      xp0(:) = xp(:,iAtom)
    
      at_norm2=SUM( at(1:3,1:3)**2 )
      IF (at_norm2.LE.Distance_Zero2) THEN
              WRITE(0,'(a)') 'You need to define lattice vectors at(1:3,i) to be&
                          & able to calculate distances'
              STOP "< Calculate_Distances >"
      END IF
      CALL Mat3Inv(at,inv_at)

      rNeigh2 = rNeigh**2
      nNeigh = 0
      dNeigh(:) = 0.d0

      DO i=1, im
         IF (i.EQ. iAtom) Cycle
         dxp(:) = xp(:,i) - xp0(:)
         ! Apply periodic boundary conditions to atom coordinates
         ! (works only ! for ~orthogonal vectors)
         dxc(:) = MatMul( inv_at(:,:), dxp(:) )
         dxc(:) = aNInt(dxc(:) )
         dxp(:) = dxp(:) - MatMul( at(:,:), dxc(:) )
         d2 = Sum( dxp(:)**2 )
         IF (d2.LE.rNeigh2) THEN
                 nNeigh = nNeigh + 1
                 dNeigh(nNeigh) = sqrt(d2)
         END IF
      END DO

      ! Sort distances in ascending order
      DO n=1, nNeigh
         nMinVec(:) = MinLoc( dNeigh(n:nNeigh) )
         nMin = nMinVec(1) + n -1
         dMin = dNeigh( nMin )
         d0 = dNeigh(n)
         dNeigh(n) = dMin
         dNeigh(nMin) = d0
      END DO

  END SUBROUTINE Calculate_Distances

END MODULE Distances_module
