PROGRAM AssociateConf

  USE associate_data
  USE associate_lecture
  USE associate_module
  USE Structure_module
  IMPLICIT NONE

  INTEGER, parameter :: verbosity_max=1
  
  CHARACTER(len=6) :: fileType
  INTEGER, EXTERNAL :: iArgc
  INTEGER, dimension(:,:), allocatable :: aux_int
  CHARACTER(len=50), dimension(1) :: aux_title
  
  ! Read name of the input file
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF

  ! Read input parameters in file 'input.dat'
  CALL Read_Associate(6)

  ! Associate atoms 
  ALLOCATE(ind(1:imm), ind0(1:imm))
  CALL Associate_Atoms(6)

  ! ==== Output structure ========================
  IF (verbosity.GE.verbosity_max) THEN
          WRITE(6,*)
          IF (outXyz) THEN
                  WRITE(6,'(2a)') 'Write output structure with Xyz format in file ', outFile
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': atom label'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': x coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': y coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 4, ': z coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 5, ': atom index'
          END IF
          IF (outCfg) THEN
                  WRITE(6,'(2a)') 'Write output structure with Cfg format in file ', outFile
                  WRITE(6,'(a,i2,a)')  '   property ', 0, ': atom index'
          END IF
          IF (outGin) WRITE(6,'(2a)') 'Write output structure with Gin format in file ', outFile
          IF (outLisa) WRITE(6,'(2a)') 'Write output structure with Lisa format in file ', outFile
          IF (outSiesta) WRITE(6,'(2a)') 'Write output structure with Siesta format in file ', outFile
          WRITE(6,*)
  END IF
  IF (outXyz.OR.outCfg.OR.outGin.OR.outLisa.OR.outSiesta) THEN

          IF (out_sort) CALL Sort_Atoms(xp, iTyp, im, ind)

          IF (outXyz) fileType='xyz'
          IF (outCfg) fileType='cfg'
          IF (outGin) fileType='gin'
          IF (outSiesta) fileType='siesta'
          IF (outLisa) fileType='lisa'

          IF (out_index) THEN
                  ALLOCATE(aux_int(1,1:imm))
                  aux_int(1,:)=ind(:)
                  aux_title(1)='atom index'
                  CALL WriteStructure(xp, iTyp, im, at, outFile, fileType, &
                        nAux_int=1, aux_int=aux_int, aux_title=aux_title)
                  DEALLOCATE(aux_int)
          ELSE
                  CALL WriteStructure(xp, iTyp, im, at, outFile, fileType)
          END IF
  END IF

                  


END PROGRAM AssociateConf
