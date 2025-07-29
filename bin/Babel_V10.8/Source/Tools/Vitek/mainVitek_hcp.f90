PROGRAM Vitek

  IMPLICIT NONE
  ! Name of the input structure file (xyz format)
  CHARACTER(len=200) :: refFile, inpFile, scriptFileScrew, scriptFileEdge, &
        outFileEdge, outFileScrew, pointFile, temp
  !INTEGER, EXTERNAL :: iArgc
  INTEGER :: im, imm, i, j, io, refUnit, inpUnit
  REAL(kind(0.d0)) :: alat, bz, R2, inv_R, Rcut, Rcut2, a0, dUz, lx, ly, x0, y0, z0, x1, y1
  REAL(kind(0.d0)), dimension(1:3,1:3) :: at
  REAL(kind(0.d0)), dimension(1:2) :: dRij, dUr
  REAL(kind(0.d0)), dimension(1:3) :: r
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp, u
  CHARACTER(len=5) :: label_temp

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter :: Distance_Zero2=Distance_Zero*Distance_Zero

  LOGICAL :: ok, lPerioZ

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
  ! Read lattice parameter
  IF (iArgc().LE.2) THEN
          WRITE(6,'(a)') 'Lattice parameter (A):'
          READ(5,*) a0
  ELSE
          CALL getArg(3,temp)
          READ(temp,*,iostat=io) a0
          IF (io.NE.0) THEN
                  WRITE(0,'(a)') 'Error when reading lattice parameter from command line arguments'
                  WRITE(0,'(2a)') '  input: ', temp
                  STOP
          END IF
  END IF

  ! Periodic boundary conditions along z
  lPerioZ=.false.
  IF (iArgc().GE.4) THEN
          CALL getArg(4,temp)
          IF ( temp(1:5).EQ."perio") lPerioZ=.true.
  END IF


  ! ---- Output files ----
  scriptFileScrew='vitek.gnu'
  outFileScrew='vitek.res'
  pointFile='vitek.points'
  scriptFileEdge='vitek_edge.gnu'
  outFileEdge='vitek_edge.res'

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


  ! Dimension of the simulation box
  x0 = MinVal(xp(1,1:im))
  y0 = MinVal(xp(2,1:im))
  x1 = MaxVal(xp(1,1:im))
  y1 = MaxVal(xp(2,1:im))
  lx = x1-x0
  ly = y1-y0

  ! ==== Periodic boundary condition along z ==
  IF (lPerioZ) THEN
          DO i=1, im
             xp(3,i) = xp(3,i) - Floor( xp(3,i)/bz )*bz
          END DO
  END IF

  ! ==== Output file ==========================

  ! ---- File containing points of the different planes ------
  OPEN(file=pointFile,unit=62,action='write',status='unknown')
  ! Atom positions corresponding to 1st (11-20) plane
  z0 = MinVal( xp(3,1:im) )
  DO i=1, im
     IF (xp(3,i).LE.z0+(0.5d0-Distance_Zero)*bz) &
             WRITE(62,'(3g14.6)') xp(1:3,i)
  END DO
  WRITE(62,*) 
  ! Atom positions corresponding to 2nd (11-20) plane
  z0 = MinVal( xp(3,1:im) ) + 0.5d0*bz
  DO i=1, im
     IF ( ( xp(3,i).GE.z0 ).AND.( xp(3,i).LE.z0+(0.5d0-Distance_Zero)*bz ) ) &
             WRITE(62,'(3g14.6)') xp(1:3,i)
  END DO
  WRITE(62,*)
  CLOSE(62)

  ! ---- Output file for screw and edge components ----------------
  OPEN(file=outFileScrew,unit=63, action="write", status="unknown")
  WRITE(63,'(a,3g14.6,a)') '#', 0.d0, 0.d0, bz, ' (Burgers vector)'
  WRITE(63,'(a,3g14.6,a)') '#', 0.d0, 0.d0, 1.d0, ' (line direction)'
  WRITE(63,'(a)') '# 1-3: 1st atom x and y coordinates'
  WRITE(63,'(a)') '# 4-6: 2nd atom x and y coordinates'
  WRITE(63,'(a)') '#   7: displacement difference uz2-uz1 along b'
  WRITE(63,'(a)') '# 8-9: atom indexes i and j'
  OPEN(file=outFileEdge,unit=64, action="write", status="unknown")
  WRITE(64,'(a,3g14.6,a)') '#', 0.d0, 0.d0, bz, ' (Burgers vector)'
  WRITE(64,'(a,3g14.6,a)') '#', 0.d0, 0.d0, 1.d0, ' (line direction)'
  WRITE(64,'(a)') '# 1-3: 1st atom x and y coordinates'
  WRITE(64,'(a)') '# 4-6: 2nd atom x and y coordinates'
  WRITE(64,'(a)') '# 7-8: displacement difference u2-u1 along x and y'
  WRITE(64,'(a)') '# 9-10: atom indexes i and j'
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
      inv_R = 1.d0/sqrt(R2)
      dUr(1:2) = u(1:2,j) - u(1:2,i)
      dUz = u(3,j) - u(3,i)
      dUz = dUz - bz*anInt(dUz/bz)
      WRITE(63,'(7g14.6,2(1x,i0))') xp(1:3,i), xp(1:3,j), dUz, i, j
      WRITE(64,'(8g14.6,2(1x,i0))') xp(1:3,i), xp(1:3,j), dUr(1:2), i, j
    END DO
  END DO
  CLOSE(63) ; CLOSE(64)

  DEALLOCATE(xp) ; DEALLOCATE(u)

  ! Gnuplot script file for screw component
  INQUIRE(file=scriptFileScrew, exist=ok)
  IF (ok) THEN
          WRITE(6,*)
          WRITE(6,'(3a)') 'Gnuplot script file ', Trim(scriptFileScrew), ' already exists'
          WRITE(6,'(a)')  'You need to delete it first, if you want a new file to be generated'
          WRITE(6,*)
  ELSE
          OPEN(file=scriptFileScrew,unit=60, action="write", status="new")
          WRITE(60,'(a)') '# Lattice parameter'
          WRITE(60,'(a,g20.12,a)') 'a0 = ', a0, '   # A'
          WRITE(60,'(a)')
          WRITE(60,'(a)') '# Burgers vector'
          WRITE(60,'(a)') 'b = a0'
          WRITE(60,'(a)')
          WRITE(60,'(a)') "# Factor for hcp structure"
          WRITE(60,'(a)') "factor=2./b"
          WRITE(60,'(a)')
          WRITE(60,'(a)') "# Threshold in A for displacement"
          WRITE(60,'(a)') "threshold = 0.1*b"
          WRITE(60,'(a)')
          WRITE(60,'(a)') "scaling(x) = abs(x)>threshold ? x*factor : 1/0"
          WRITE(60,'(a)')
          WRITE(60,'(a)') "set pointsize 1"
          WRITE(60,'(a)') "set border 0"
          WRITE(60,'(a)') "unset ytics"
          WRITE(60,'(a)') "unset xtics"
          WRITE(60,'(a,g14.6,a)') "lx = ", lx, "*1.01"
          WRITE(60,'(a,g14.6,a)') "ly = ", ly, "*1.01"
          WRITE(60,'(a,g14.6,a)') "x0 = ", x0, "-0.005*lx"
          WRITE(60,'(a,g14.6,a)') "y0 = ", y0, "-0.005*ly"
          WRITE(60,'(a)') "set xrange [x0:x0+lx]"
          WRITE(60,'(a)') "set yrange [y0:y0+ly]"
          WRITE(60,'(a)') "set size ratio -1"
          WRITE(60,*)
          WRITE(60,'(a)') "unset key"
          WRITE(60,'(3a)') "plot '", Trim(pointFile), "' u 1:2 every :2::0 w p pt 31 lc rgb 'black',\"
          WRITE(60,'(3a)') "     '", Trim(pointFile), "' u 1:2 every :2::1 w p pt 65 lc rgb 'black',\"
          WRITE(60,'(3a)') "     '", Trim(outFileScrew), "' u &
                &(0.5*($1+$4)-0.5*factor*$7*($4-$1)):(0.5*($2+$5)-0.5*factor*$7*($5-$2)):(scaling($7)*($4-$1)):(scaling($7)*($5-$2)) &
                &w vectors filled head lt 1"
          WRITE(60,*)
          WRITE(60,'(a)') "set output 'vitek.eps'"
          WRITE(60,'(a)') "set term post eps enh color solid 'Times-Roman' 20"
          WRITE(60,'(a)') "replot"
          WRITE(60,'(a)') "set output"
          WRITE(60,'(a)') "set term pop"
          WRITE(60,*)
          CLOSE(60)
  END IF

  ! Gnuplot script file for Edge component
  INQUIRE(file=scriptFileEdge, exist=ok)
  IF (ok) THEN
          WRITE(6,*)
          WRITE(6,'(3a)') 'Gnuplot script file ', Trim(scriptFileEdge), ' already exists'
          WRITE(6,'(a)')  'You need to delete it first, if you want a new file to be generated'
          WRITE(6,*)
  ELSE
          OPEN(file=scriptFileEdge,unit=61, action="write", status="new")
          WRITE(61,'(a)') "# Multiplying factor"
          WRITE(61,'(a)') "factor=1."
          WRITE(61,'(a)')
          WRITE(61,'(a)') "# Threshold in A for displacement"
          WRITE(61,'(a)') "threshold = 0."
          WRITE(61,'(a)')
          WRITE(61,'(a)') "xScaling(x,y) = sqrt(x**2+y**2)>threshold ? x*factor : 1/0"
          WRITE(61,'(a)') "yScaling(x,y) = sqrt(x**2+y**2)>threshold ? y*factor : 1/0"
          WRITE(61,'(a)')
          WRITE(61,'(a)') "set pointsize 1"
          WRITE(61,'(a)') "set border 0"
          WRITE(61,'(a)') "unset ytics"
          WRITE(61,'(a)') "unset xtics"
          WRITE(61,'(a,g14.6,a)') "lx = ", lx, "*1.01"
          WRITE(61,'(a,g14.6,a)') "ly = ", ly, "*1.01"
          WRITE(61,'(a,g14.6,a)') "x0 = ", x0, "-0.005*lx"
          WRITE(61,'(a,g14.6,a)') "y0 = ", y0, "-0.005*ly"
          WRITE(61,'(a)') "set xrange [x0:x0+lx]"
          WRITE(61,'(a)') "set yrange [y0:y0+ly]"
          WRITE(61,'(a)') "set size ratio -1"
          WRITE(61,*)
          WRITE(61,'(a)') "unset key"
          WRITE(61,'(3a)') "plot '", Trim(pointFile), "' u 1:2 every :3::0 w p pt 31 lc rgb 'black',\"
          WRITE(61,'(3a)') "     '", Trim(pointFile), "' u 1:2 every :3::1 w p pt 65 lc rgb 'black',\"
          WRITE(61,'(3a)') "     '", Trim(outFileEdge), "' u &
                &(0.5*($1+$4)-0.5*factor*$7):(0.5*($2+$5)-0.5*factor*$8):&
                &(xScaling($7,$8)):(yScaling($7,$8)) &
                &w vectors filled head lt 1"
          WRITE(61,*)
          WRITE(61,'(a)') "set output 'vitek_edge.eps'"
          WRITE(61,'(a)') "set term post eps enh solid 'Times-Roman' 20"
          WRITE(61,'(a)') "replot"
          WRITE(61,'(a)') "set output"
          WRITE(61,'(a)') "set term pop"
          WRITE(61,*)
          CLOSE(61)
  END IF

END PROGRAM Vitek
