PROGRAM Vitek

  IMPLICIT NONE
  ! Name of the input structure file (xyz format)
  CHARACTER(len=200) :: refFile, inpFile, outFileScrew, outFileEdge, temp
  INTEGER, EXTERNAL :: iArgc
  INTEGER :: im, imm, i, j, io, refUnit, inpUnit, outUnitScrew, outUnitEdge
  REAL(kind(0.d0)) :: alat, bz, R2, inv_R, Rcut, Rcut2, a0, dUz, lx, ly, x0, y0, z0, x1, y1
  REAL(kind(0.d0)), dimension(1:3,1:3) :: at
  REAL(kind(0.d0)), dimension(1:2) :: dRij, xm, dUr
  REAL(kind(0.d0)), dimension(1:3) :: r
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
  ! Read name of the output files
  IF (iArgc().LE.2) THEN
          WRITE(6,'(a)') 'Name of the output file for screw component (gnuplot script)'
          READ(5,*) outFileScrew
  ELSE
          CALL getArg(3,outFileScrew)
  END IF
  IF (iArgc().LE.3) THEN
          WRITE(6,'(a)') 'Name of the output file for edge component (gnuplot script)'
          READ(5,*) outFileEdge
  ELSE
          CALL getArg(4,outFileEdge)
  END IF
  ! Read lattice parameter
  IF (iArgc().LE.4) THEN
          WRITE(6,'(a)') 'Lattice parameter (A):'
          READ(5,*) a0
  ELSE
          CALL getArg(5,temp)
          READ(temp,*,iostat=io) a0
          IF (io.NE.0) THEN
                  WRITE(0,'(a)') 'Error when reading lattice parameter from command line arguments'
                  WRITE(0,'(2a)') '  input: ', temp
                  STOP
          END IF
  END IF
  ! Input parameters for screw dislocation in Fe
  bz = a0*sqrt(3.d0)*0.5d0
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
  outUnitScrew=60
  outUnitEdge=61
  OPEN(file=outFileScrew,unit=outUnitScrew, action="write", status="unknown")
  OPEN(file=outFileEdge,unit=outUnitEdge, action="write", status="unknown")

  ! Screw component
  WRITE(outUnitScrew,'(a)') "# Factor for bcc structure"
  WRITE(outUnitScrew,'(a)') "factor=sqrt(6.)/3./(sqrt(3.)/6)"
  WRITE(outUnitScrew,'(a)')
  WRITE(outUnitScrew,'(a)') "# Threshold in A for displacement"
  WRITE(outUnitScrew,'(a)') "threshold = 0.05"
  WRITE(outUnitScrew,'(a)')
  WRITE(outUnitScrew,'(a)') "set pointsize 1"
  WRITE(outUnitScrew,'(a)') "set border 0"
  WRITE(outUnitScrew,'(a)') "unset ytics"
  WRITE(outUnitScrew,'(a)') "unset xtics"
  WRITE(outUnitScrew,'(a,g14.6,a)') "lx = ", lx, "*1.01"
  WRITE(outUnitScrew,'(a,g14.6,a)') "ly = ", ly, "*1.01"
  WRITE(outUnitScrew,'(a,g14.6,a)') "x0 = ", x0, "-0.005*lx"
  WRITE(outUnitScrew,'(a,g14.6,a)') "y0 = ", y0, "-0.005*ly"
  WRITE(outUnitScrew,'(a)') "set xrange [x0:x0+lx]"
  WRITE(outUnitScrew,'(a)') "set yrange [y0:y0+ly]"
  WRITE(outUnitScrew,'(a)') "set size ratio ly/lx"
  WRITE(outUnitScrew,*)
  WRITE(outUnitScrew,'(a)') "# set output 'vitek_screw.jpg'"
  WRITE(outUnitScrew,'(a)') "# set term jpeg large"
  WRITE(outUnitScrew,'(a)') "set output 'vitek_screw.eps'"
  WRITE(outUnitScrew,'(a)') "set term post eps enh solid 'Times-Roman' 20"
  WRITE(outUnitScrew,*)
  WRITE(outUnitScrew,'(a)') "unset key"
  WRITE(outUnitScrew,'(a)') "plot '-' u 1:2 w p pt 31 lc rgb 'black',\"
  WRITE(outUnitScrew,'(a)') "     '-' u 1:2 w p pt 65 lc rgb 'black',\"
  WRITE(outUnitScrew,'(a)') "     '-' u 1:2 w p pt 31 lc rgb 'grey',\"
  WRITE(outUnitScrew,'(a)') "'-' u ($1-0.5*factor*$3):($2-0.5*factor*$4):&
        &(sqrt($3**2+$4**2)>threshold?factor*$3:1/0):(sqrt($3**2+$4**2)>threshold?factor*$4:1/0) &
        &w vectors filled head lt 2"

  ! Edge component
  WRITE(outUnitEdge,'(a)') "# Multiplying factor"
  WRITE(outUnitEdge,'(a)') "factor=20."
  WRITE(outUnitEdge,'(a)')
  WRITE(outUnitEdge,'(a)') "# Threshold in A for displacement"
  WRITE(outUnitEdge,'(a)') "threshold = 0."
  WRITE(outUnitEdge,'(a)')
  WRITE(outUnitEdge,'(a)') "set pointsize 1"
  WRITE(outUnitEdge,'(a)') "set border 0"
  WRITE(outUnitEdge,'(a)') "unset ytics"
  WRITE(outUnitEdge,'(a)') "unset xtics"
  WRITE(outUnitEdge,'(a,g14.6,a)') "lx = ", lx, "*1.01"
  WRITE(outUnitEdge,'(a,g14.6,a)') "ly = ", ly, "*1.01"
  WRITE(outUnitEdge,'(a,g14.6,a)') "x0 = ", x0, "-0.005*lx"
  WRITE(outUnitEdge,'(a,g14.6,a)') "y0 = ", y0, "-0.005*ly"
  WRITE(outUnitEdge,'(a)') "set xrange [x0:x0+lx]"
  WRITE(outUnitEdge,'(a)') "set yrange [y0:y0+ly]"
  WRITE(outUnitEdge,'(a)') "set size ratio ly/lx"
  WRITE(outUnitEdge,*)
  WRITE(outUnitEdge,'(a)') "# set output 'vitek_edge.jpg'"
  WRITE(outUnitEdge,'(a)') "# set term jpeg large"
  WRITE(outUnitEdge,'(a)') "set output 'vitek_edge.eps'"
  WRITE(outUnitEdge,'(a)') "set term post eps enh solid 'Times-Roman' 20"
  WRITE(outUnitEdge,*)
  WRITE(outUnitEdge,'(a)') "unset key"
  WRITE(outUnitEdge,'(a)') "plot '-' u 1:2 w p pt 31 lc rgb 'black',\"
  WRITE(outUnitEdge,'(a)') "     '-' u 1:2 w p pt 65 lc rgb 'black',\"
  WRITE(outUnitEdge,'(a)') "     '-' u 1:2 w p pt 31 lc rgb 'grey',\"
  WRITE(outUnitEdge,'(a)') "'-' u ($1-0.5*factor*$3):($2-0.5*factor*$4):&
        &(sqrt($3**2+$4**2)>threshold?factor*$3:1/0):(sqrt($3**2+$4**2)>threshold?factor*$4:1/0) &
        &w vectors filled heads lt 2"

  z0 = Sum( Modulo(xp(3,1:im),bz) )/dble(im)-bz/6.
  ! Atom positions corresponding to 1st (111) plane
  DO i=1, im
     IF (3.d0*Modulo(xp(3,i)-z0,bz).LT.bz) THEN
             WRITE(outUnitScrew,'(3g14.6)') xp(1:3,i)
             WRITE(outUnitEdge,'(3g14.6)') xp(1:3,i)
     END IF
  END DO
  WRITE(outUnitScrew,'(a)') "e"
  WRITE(outUnitEdge,'(a)') "e"
  ! Atom positions corresponding to 2nd (111) plane
  DO i=1, im
     IF ( (3.d0*Modulo(xp(3,i)-z0,bz).GE.bz).AND. (3.d0*Modulo(xp(3,i)-z0,bz).LT.2.d0*bz) ) THEN
             WRITE(outUnitScrew,'(3g14.6)') xp(1:3,i)
             WRITE(outUnitEdge,'(3g14.6)') xp(1:3,i)
     END IF
  END DO
  WRITE(outUnitScrew,'(a)') "e"
  WRITE(outUnitEdge,'(a)') "e"
  ! Atom positions corresponding to 3rd (111) plane
  DO i=1, im
     IF (3.d0*Modulo(xp(3,i)-z0,bz).GE.2.d0*bz) THEN
             WRITE(outUnitScrew,'(3g14.6)') xp(1:3,i)
             WRITE(outUnitEdge,'(3g14.6)') xp(1:3,i)
     END IF
  END DO
  WRITE(outUnitScrew,'(a)') "e"
  WRITE(outUnitEdge,'(a)') "e"

  ! Arrow definitions
  z0 = MinVal( xp(3,1:im)) + 2.d0*bz
  DO i=1, im-1
    IF (xp(3,i).GT.z0) Cycle
    DO j=i+1, im
      IF (xp(3,j).GT.z0) Cycle
      dRij(:) = xp(1:2,j) - xp(1:2,i) 
      R2 = dRij(1)**2+dRij(2)**2
      IF (R2.GT.Rcut2) Cycle
      IF (R2.LE.distance_zero2) Cycle
      inv_R = 1.d0/sqrt(R2)
      dUr(1:2) = u(1:2,j) - u(1:2,i)
      dUz = u(3,j) - u(3,i)
      dUz = dUz - bz*anInt(dUz/bz)
      xm(1:2) = 0.5d0*( xp(1:2,j) +xp(1:2,i) )
      WRITE(outUnitScrew,'(4g14.6)') xm(1:2), dUz*inv_R*dRij(1:2)
      WRITE(outUnitEdge,'(4g14.6)') xm(1:2), dUr(1:2)
    END DO
  END DO
  WRITE(outUnitScrew,'(a)') "e"
  WRITE(outUnitScrew,*)
  WRITE(outUnitScrew,'(a)') "set output"
  WRITE(outUnitScrew,'(a)') "set term x11"
  WRITE(outUnitEdge,'(a)') "e"
  WRITE(outUnitEdge,*)


  IF (outUnitScrew.NE.6) CLOSE(outUnitScrew)
  IF (outUnitEdge.NE.6) CLOSE(outUnitEdge)
  DEALLOCATE(xp) ; DEALLOCATE(u)

END PROGRAM Vitek
