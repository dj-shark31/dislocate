PROGRAM VitekExtract

  IMPLICIT NONE
  ! Name of the input structure file (xyz format)
  CHARACTER(len=200) :: refFile, inpFile, outFile, temp
  INTEGER, EXTERNAL :: iArgc
  INTEGER :: im, imm, i, j, io, refUnit, inpUnit, outUnit
  REAL(kind(0.d0)) :: alat, bz, R2, Rcut, Rcut2, a0, z0, dUz
  REAL(kind(0.d0)), dimension(1:3,1:3) :: at
  REAL(kind(0.d0)), dimension(1:3) :: r
  REAL(kind(0.d0)), dimension(1:2) :: dRij
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp, u
  CHARACTER(len=5) :: label_temp

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter :: Distance_Zero2=Distance_Zero*Distance_Zero


  ! ==== Command line arguments ================
  ! Read name of the reference file containing atom positions
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the reference file with atom positions (xyz format)'
          READ(5,*) refFile
  ELSE
          CALL getArg(1,refFile)
  END IF
  ! Read name of the input file containing atom displacement
  IF (iArgc().LE.1) THEN
          WRITE(6,'(a)') 'Name of the input file with atom displacements (xyz format)'
          READ(5,*) inpFile
  ELSE
          CALL getArg(2,inpFile)
  END IF
  ! Read name of the output file
  IF (iArgc().LE.2) THEN
          WRITE(6,'(a)') 'Name of the output file containing displacement differences'
          READ(5,*) outFile
  ELSE
          CALL getArg(3,outFile)
  END IF
  ! Read lattice parameter
  IF (iArgc().LE.3) THEN
          WRITE(6,'(a)') 'Lattice parameter (A):'
          READ(5,*) a0
  ELSE
          CALL getArg(4,temp)
          READ(temp,*,iostat=io) a0
          IF (io.NE.0) THEN
                  WRITE(0,'(a)') 'Error when reading lattice parameter from command line arguments'
                  WRITE(0,'(2a)') '  input: ', temp
                  STOP
          END IF
  END IF
  ! Input parameters for a screw dislocation in hcp lattice
  bz = a0
  Rcut = 0.5d0*a0*( 1.d0 + sqrt(3.d0) )  ! Half distance between 1st and 2nd n.n.
  Rcut2 = Rcut**2

  ! ===== Reference structure containing atom positions =====================
  refUnit=50
  OPEN(file=refFile,unit=refUnit, action="read", status="old")

  READ(refUnit,*) im ! Number of atoms
  READ(refUnit,*)    ! Line for title

  imm=im        ! Maximal number of atoms for table allocations
  ALLOCATE(xp(1:3,1:imm)) ; ALLOCATE(u(1:3,1:imm))

  ! Read atom coordinates 
  DO i=1, im
     READ(refUnit,*) label_temp, xp(1:3,i)
  END DO

  ! Read periodicity vectors if defined
  at(1:3,1:3) = 0.d0
  READ(refUnit,*,iostat=io)
  IF (io.EQ.0) THEN
          DO i=1, 3
             READ(refUnit,*,iostat=io) at(1:3,i)
             IF (io.NE.0) Exit
          END DO
  END IF
  IF (io.NE.0) at(1:3,1:3) = 0.d0

  ! Read scaling factor if defined
  alat=1.d0
  IF (io.EQ.0) READ(refUnit,*,iostat=io)
  IF (io.EQ.0) READ(refUnit,*,iostat=io) alat
  IF (io.NE.0) alat=1.D0

  IF (refUnit.NE.5) CLOSE(refUnit) 

  ! Rescale atom coordinates and displacements
  xp(:,1:im) = alat*xp(:,1:im)

  ! ===== Input structure containing atom displacement =====================
  inpUnit=51
  OPEN(file=inpFile,unit=inpUnit, action="read", status="old")

  READ(inpUnit,*) im ! Number of atoms
  READ(inpUnit,*)    ! Line for title

  IF (im.NE.imm) THEN
          WRITE(0,'(a,i0)') 'Number of atoms read in reference structure file: ', imm
          WRITE(0,'(a,i0)') 'Number of atoms read in input structure file: ', im
          STOP '< Vitek >'
  END IF

  ! Read atom coordinates and displacement vector
  DO i=1, im
     READ(inpUnit,*,iostat=io) label_temp, r(1:3), u(1:3,i)
     IF (io.NE.0) THEN
             WRITE(0,'(2a)') 'WARNING: Error when reading atom displacement in file ', Trim(inpFile)
             WRITE(0,'(a)') '   => null displacement is assumed for all atoms'
             u(:,:)=0.d0
             EXIT
     END IF
  END DO

  ! Read periodicity vectors if defined
  at(1:3,1:3) = 0.d0
  READ(inpUnit,*,iostat=io)
  IF (io.EQ.0) THEN
          DO i=1, 3
             READ(inpUnit,*,iostat=io) at(1:3,i)
             IF (io.NE.0) Exit
          END DO
  END IF
  IF (io.NE.0) at(1:3,1:3) = 0.d0

  ! Read scaling factor if defined
  alat=1.d0
  IF (io.EQ.0) READ(inpUnit,*,iostat=io)
  IF (io.EQ.0) READ(inpUnit,*,iostat=io) alat
  IF (io.NE.0) alat=1.D0

  IF (inpUnit.NE.5) CLOSE(inpUnit) 

  ! Rescale atom displacements
  u(:,1:im) = alat*u(:,1:im)

  ! ==== Output file ==========================
  outUnit=62
  OPEN(file=outFile,unit=outUnit, action="write", status="unknown")
  WRITE(outUnit,'(a,3g14.6,a)') '#', 0.d0, 0.d0, bz, ' (Burgers vector)'
  WRITE(outUnit,'(a,3g14.6,a)') '#', 0.d0, 0.d0, 1.d0, ' (line direction)'
  WRITE(outUnit,'(a)') '# 1-3: 1st atom x and y coordinates'
  WRITE(outUnit,'(a)') '# 4-6: 2nd atom x and y coordinates'
  WRITE(outUnit,'(a)') '#   7: displacement difference uz2-uz1 along b'

  ! Arrow definitions
  z0 = MinVal( xp(3,1:im)) + (1.d0-Distance_Zero)*bz
  DO i=1, im-1
    IF (xp(3,i).GT.z0) Cycle
    DO j=i+1, im
      IF (xp(3,j).GT.z0) Cycle
      dRij(:) = xp(1:2,j) - xp(1:2,i) 
      R2 = dRij(1)**2+dRij(2)**2
      IF (R2.GT.Rcut2) Cycle
      IF (R2.LE.distance_zero2) Cycle
      dUz = u(3,j) - u(3,i)
      dUz = dUz - bz*anInt(dUz/bz)
      WRITE(outUnit,'(7g14.6)') xp(1:3,i), xp(1:3,j), dUz
    END DO
  END DO

  IF (outUnit.NE.6) CLOSE(outUnit)
  DEALLOCATE(xp) ; DEALLOCATE(u)

END PROGRAM VitekExtract
