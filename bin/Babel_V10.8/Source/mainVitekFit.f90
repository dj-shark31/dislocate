PROGRAM VitekFit

  USE babel_data
  USE Disloc_elasticity_ani
  USE lineCouple_elasticity_ani
  USE lineForce_elasticity_ani
  USE DDipoleModule
  USE babel_lecture
  USE slabModule
  USE loopModule
  USE fitModule
  USE fitVitekModule
  USE fitConstraintModule
  IMPLICIT NONE
 
  INTEGER, parameter :: verbosity_max=1
  !INTEGER, EXTERNAL :: iArgc

  ! Name of the file containing Vitek map
  CHARACTER(len=100) :: Vitek_file

  program_name="vitekFit"

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
  ! Read name of the file containing Vitek map
  IF (iArgc().LE.1) THEN
          WRITE(6,'(a)') 'Name of the Vitek file'
          READ(5,*) Vitek_file
  ELSE
          CALL getArg(2,Vitek_file)
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
  IF (duplicate)  STOP 'Option not supported'
  IF (rotate)     STOP 'Option not supported'
  IF (translate)  STOP 'Option not supported'
  IF (slab)       STOP 'Option not supported'
  IF (remove_cut) STOP 'Option not supported'
  IF (l_loop)     STOP 'Option not supported'
  IF (lineForce)  STOP 'Option not supported'
  IF (l_lineCouple) STOP 'Option not supported'

  ! Read definition of Vitek map  
  OPEN(file=Vitek_file,unit=50,status='old',action='read')
  CALL ReadVitek(50, threshold)
  CLOSE(50)

  ! Fit dislocation positions to Vitek map
  CALL FitVitek(6)

  ! Write results
  IF (l_dislo) CALL PrintDislo(6)
  IF (l_DDipole) THEN
         CALL PrintDDipole(6)
         OPEN(file="vitekPara.res", unit=60, action='write')
         CALL WriteDDipole(60) 
         CLOSE(60)
 END IF


END PROGRAM VitekFit
