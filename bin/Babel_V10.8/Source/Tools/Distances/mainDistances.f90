PROGRAM Distances

  USE structure_module
  USE distances_module

  IMPLICIT NONE
  
  INTEGER, EXTERNAL :: iArgc
  REAL(kind(0.d0)), dimension(:), allocatable :: dNeigh
  INTEGER :: nNeigh, n

  ! Read name of the input file
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF

  ! Read input file
  CALL Read_distances(6)

  ! Duplicate unit cell
  IF (duplicate) THEN
          CALL duplicate_cell(lat, xp, iTyp, im, at, 6, verbose=verbosity.GE.3)
          duplicate=.FALSE.
          Lat(1:3)=1
  END IF

  ! Calculate distance with neighbours
  ALLOCATE(dNeigh(1:imm))
  CALL Calculate_Distances(dNeigh, nNeigh)

  WRITE(6,*)
  WRITE(6,'(a)') "####################################"
  WRITE(6,*)
  WRITE(6,'(a,i0)') "# Distances for atom ", iAtom
  WRITE(6,'(a,g20.12)') "#  cutoff distance, rNeigh = ", rNeigh
  DO n=1, nNeigh
     WRITE(6,'(g20.12)') dNeigh(n)
  END DO



END PROGRAM Distances
