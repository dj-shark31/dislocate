PROGRAM Vitek

  IMPLICIT NONE
  ! Name of the input structure file (xyz format)
  CHARACTER(len=200) :: refFile, inpFile, outFile
  INTEGER, EXTERNAL :: iArgc
  INTEGER :: im, imm, i, j, io, refUnit, inpUnit, outUnit
  REAL(kind(0.d0)) :: alat, b, R2, inv_R, Rcut, Rcut2, a0, dU, lx, ly, x0, y0, z0, x1, y1
  REAL(kind(0.d0)), dimension(1:3,1:3) :: at
  REAL(kind(0.d0)), dimension(1:2) :: dRij, xm
  REAL(kind(0.d0)), dimension(1:3) :: r, bVector
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp, u
  CHARACTER(len=5) :: label_temp

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter :: Distance_Zero2=Distance_Zero*Distance_Zero


  ! ==== Command line arguments ================
  ! Read name of the reference file containing atom positions
  WRITE(6,'(a)') 'Name of the reference file with atom positions (xyz format)'
  READ(5,*) refFile
  ! Read name of the input file containing atom displacement
  WRITE(6,'(a)') 'Name of the input file with atom displacements (xyz format)'
  READ(5,*) inpFile
  ! Read name of the output files
  WRITE(6,'(a)') 'Name of the output file (gnuplot script)'
  READ(5,*) outFile
  ! Read lattice parameter
  WRITE(6,'(a)') 'Lattice parameter (A):'
  READ(5,*) a0
  ! Read Burgers vector direction
  WRITE(6,'(a)') 'Direction for projecting displacement (Burgers vector direction):'
  READ(5,*) bVector(:)
  bVector(:) = bVector(:)/Sqrt( Sum( bVector(:)**2 ) )
  ! Input parameters for dislocation in bcc crystal
  b = a0*sqrt(3.d0)*0.5d0
  Rcut = 0.5d0*a0*( sqrt(6.d0)/3.d0 + sqrt(2.d0) )  ! Half distance between 1st and 2nd n.n.
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


  ! Dimension of the simulation box
  x0 = MinVal(xp(1,1:im))
  y0 = MinVal(xp(2,1:im))
  x1 = MaxVal(xp(1,1:im))
  y1 = MaxVal(xp(2,1:im))
  lx = x1-x0
  ly = y1-y0

  ! ==== Output file ==========================
  outUnit=60
  OPEN(file=outFile,unit=outUnit, action="write", status="unknown")

  ! Component along Burgers vector
  WRITE(outUnit,'(a)') "# Factor for bcc structure"
  WRITE(outUnit,'(a)') "factor=sqrt(6.)/3./(sqrt(3.)/6)"
  WRITE(outUnit,'(a)')
  WRITE(outUnit,'(a)') "# Threshold in A for displacement"
  WRITE(outUnit,'(a)') "threshold = 0.05"
  WRITE(outUnit,'(a)')
  WRITE(outUnit,'(a)') "set pointsize 1"
  WRITE(outUnit,'(a)') "set border 0"
  WRITE(outUnit,'(a)') "unset ytics"
  WRITE(outUnit,'(a)') "unset xtics"
  WRITE(outUnit,'(a,g14.6,a)') "lx = ", lx, "*1.01"
  WRITE(outUnit,'(a,g14.6,a)') "ly = ", ly, "*1.01"
  WRITE(outUnit,'(a,g14.6,a)') "x0 = ", x0, "-0.005*lx"
  WRITE(outUnit,'(a,g14.6,a)') "y0 = ", y0, "-0.005*ly"
  WRITE(outUnit,'(a)') "set xrange [x0:x0+lx]"
  WRITE(outUnit,'(a)') "set yrange [y0:y0+ly]"
  WRITE(outUnit,'(a)') "set size ratio ly/lx"
  WRITE(outUnit,*)
  WRITE(outUnit,'(a)') "# set output 'vitek.jpg'"
  WRITE(outUnit,'(a)') "# set term jpeg large"
  WRITE(outUnit,'(a)') "set output 'vitek.eps'"
  WRITE(outUnit,'(a)') "set term post eps enh solid 'Times-Roman' 20"
  WRITE(outUnit,*)
  WRITE(outUnit,'(a)') "unset key"
  WRITE(outUnit,'(a)') "plot '-' u 1:2 w p pt 31 lc rgb 'black',\"
  WRITE(outUnit,'(a)') "     '-' u 1:2 w p pt 65 lc rgb 'black',\"
  WRITE(outUnit,'(a)') "     '-' u 1:2 w p pt 31 lc rgb 'grey',\"
  WRITE(outUnit,'(a)') "'-' u ($1-0.5*factor*$3):($2-0.5*factor*$4):&
        &(sqrt($3**2+$4**2)>threshold?factor*$3:1/0):(sqrt($3**2+$4**2)>threshold?factor*$4:1/0) &
        &w vectors filled head lt 2"


  z0 = Sum( Modulo(xp(3,1:im),b) )/dble(im)-b/6.
  ! Atom positions corresponding to 1st (111) plane
  DO i=1, im
     IF (3.d0*Modulo(xp(3,i)-z0,b).LT.b) THEN
             WRITE(outUnit,'(3g14.6)') xp(1:3,i)
     END IF
  END DO
  WRITE(outUnit,'(a)') "e"
  ! Atom positions corresponding to 2nd (111) plane
  DO i=1, im
     IF ( (3.d0*Modulo(xp(3,i)-z0,b).GE.b).AND. (3.d0*Modulo(xp(3,i)-z0,b).LT.2.d0*b) ) THEN
             WRITE(outUnit,'(3g14.6)') xp(1:3,i)
     END IF
  END DO
  WRITE(outUnit,'(a)') "e"
  ! Atom positions corresponding to 3rd (111) plane
  DO i=1, im
     IF (3.d0*Modulo(xp(3,i)-z0,b).GE.2.d0*b) THEN
             WRITE(outUnit,'(3g14.6)') xp(1:3,i)
     END IF
  END DO
  WRITE(outUnit,'(a)') "e"

  ! Arrow definitions
  z0 = MinVal( xp(3,1:im)) + 2.d0*b
  DO i=1, im-1
    IF (xp(3,i).GT.z0) Cycle
    DO j=i+1, im
      IF (xp(3,j).GT.z0) Cycle
      dRij(:) = xp(1:2,j) - xp(1:2,i) 
      R2 = dRij(1)**2+dRij(2)**2
      IF (R2.GT.Rcut2) Cycle
      IF (R2.LE.distance_zero2) Cycle
      inv_R = 1.d0/sqrt(R2)
      dU = Sum( ( u(:,j) - u(:,i) )*bVector(:) )
      dU = dU - b*anInt(dU/b)
      xm(1:2) = 0.5d0*( xp(1:2,j) +xp(1:2,i) )
      WRITE(outUnit,'(4g14.6)') xm(1:2), dU*inv_R*dRij(1:2)
    END DO
  END DO
  WRITE(outUnit,'(a)') "e"
  WRITE(outUnit,*)

  IF (outUnit.NE.6) CLOSE(outUnit)
  DEALLOCATE(xp) ; DEALLOCATE(u)

END PROGRAM Vitek
