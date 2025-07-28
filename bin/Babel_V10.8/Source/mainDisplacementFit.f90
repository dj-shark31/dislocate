PROGRAM DisplacementFit

  USE babel_data
  USE Disloc_elasticity_ani
  USE lineCouple_elasticity_ani
  USE lineForce_elasticity_ani
  USE DDipoleModule
  USE babel_lecture
  USE slabModule
  USE loopModule
  USE fitModule
  USE fitDisplacementModule
  USE fitConstraintModule
  IMPLICIT NONE
 
  INTEGER, parameter :: verbosity_max=1
  !INTEGER, EXTERNAL :: iArgc

  ! Name of the file containing displacement map
  CHARACTER(len=100) :: displacement_file

  program_name="displacementFit"

  ! Initialize quantities controlling the fit
  CALL InitFit()
  CALL InitConstraint()

  ! Read name of the input file
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF
  ! Read name of the file containing displacement map
  IF (iArgc().LE.1) THEN
          WRITE(6,'(a)') 'Name of the displacement file'
          READ(5,*) displacement_file
  ELSE
          CALL getArg(2,displacement_file)
  END IF

  ! Read input parameters in file 'input.dat'
  CALL Read_Babel(6)

  ! Print fitting parameters
  CALL PrintFit(6)
  CALL PrintConstraint(6)

  ! Reinitialize atom coordinates, keeping only periodicity vectors
  im=0 ; im0=0 ; im1=0
  IF (Allocated(xp)) DEALLOCATE(xp)
  IF (Allocated(xp0)) DEALLOCATE(xp0)
  IF (Allocated(xp1)) DEALLOCATE(xp1)
  IF (Allocated(iTyp)) DEALLOCATE(iTyp)
  IF (Allocated(iTyp0)) DEALLOCATE(iTyp0)
  IF (Allocated(iTyp1)) DEALLOCATE(iTyp1)
  IF (Allocated(keep)) DEALLOCATE(keep)

  ! Babel options which are not supported
  IF (duplicate) STOP 'Option not supported'
  IF (rotate) STOP 'Option not supported'
  IF (translate) STOP 'Option not supported'
  IF (slab) STOP 'Option not supported'
  IF (remove_cut) STOP 'Option not supported'
  IF (l_loop) STOP 'Option not supported'
  IF (lineForce) STOP 'Option not supported'
  IF (l_lineCouple) STOP 'Option not supported'

  ! Read definition of displacement map  
  OPEN(file=displacement_file,unit=50,status='old',action='read')
  CALL ReadDisplacement(50)
  CLOSE(50)

  ! Fit dislocation positions to Displacement map
  CALL FitDisplacement(6)

  ! Write results
  IF (l_dislo) CALL PrintDislo(6)
  IF (l_DDipole) THEN
         CALL PrintDDipole(6)
         OPEN(file="displacementPara.res", unit=60, action='write')
         CALL WriteDDipole(60) 
         CLOSE(60)
 END IF


END PROGRAM DisplacementFit
