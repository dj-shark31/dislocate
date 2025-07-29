PROGRAM PatternInit

  USE gradElasticDisplacementModule

  IMPLICIT NONE
 
  INTEGER, parameter :: verbosity_max=1
  INTEGER :: verbosity
  LOGICAL :: debug

  REAL(kind(0.d0)) :: alat     ! Lattice parameter (A)
  REAL(kind(0.d0)) :: coa      ! c/a ratio for hcp
  CHARACTER(len=3) :: crystal   ! hcp/fcc/bcc
  LOGICAL :: addHcpBasalFault   ! with patterns corresponding to basal SF (hcp)
  LOGICAL :: addHcpPrismFault   ! with patterns corresponding to prism SF (hcp)
  LOGICAL :: addHcpPy1Fault     ! with patterns corresponding to pyramidal 1 SF (hcp)
  LOGICAL :: addHcpPy1Twin      ! with patterns corresponding to 1st order pyramidal twins (hcp)
  LOGICAL :: addHcpPy1TB        ! with patterns corresponding to 1st order pyramidal twin boundaries (hcp)
  LOGICAL :: rotate             ! rotation
  REAL(kind(0.d0)), dimension(3,3) :: rot
  !INTEGER, EXTERNAL :: iArgc

  ! Name of the input file defining the dislocation
  CHARACTER(len=100) :: input_file

  ! Read name of the input file
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF

  ! Read input parameters in file 'input.dat'
  CALL read_inputPattern(6)

  SELECT CASE(Trim(crystal)) 
  CASE('hcp') 
          CALL InitHcpPattern(alat,coa)
          IF (addHcpBasalFault) CALL HcpPattern_AddBasalFault(alat,coa)
          IF (addHcpPy1Fault)   CALL HcpPattern_AddPy1Fault(alat,coa)
          IF (addHcpPrismFault) CALL HcpPattern_AddPrismFault(alat,coa)
          IF (addHcpPy1Twin)    CALL HcpPattern_AddPy1Twin(alat,coa)
          IF (addHcpPy1TB)      CALL HcpPattern_AddPy1TB(alat,coa)
  CASE('bcc') 
          CALL InitBccPattern(alat)
  CASE DEFAULT
          WRITE(6,'(2a)') 'Unknown crystal type: crystal = ', Trim(crystal)
  END SELECT

  ! Rotation
  IF (rotate) THEN  
          CALL rotatePattern(rot) 
          rotate=.FALSE.
  END IF

  ! Print pattern definition
  IF (verbosity.GE.verbosity_max) CALL PrintAllPatterns(6)

  ! Write pattern definition in a file
  IF (patternFile.NE.'') THEN
          OPEN(file=patternFile,unit=60,action='write')
          CALL WritePattern(60)
          CLOSE(60)
          IF (verbosity.GE.verbosity_max) THEN
                  WRITE(6,*)
                  WRITE(6,'(2a)') 'Pattern definition written in file ', Trim(patternFile)
                  WRITE(6,*)
          END IF
  END IF

CONTAINS

  SUBROUTINE read_inputPattern(out)
  
      USE GradElasticDisplacementModule, ONLY : patternFile
      IMPLICIT NONE
      INTEGER, intent(in) :: out
  
      INTEGER, parameter :: verbosity_max=1

      NAMELIST /input/ crystal, alat, coa, &
          addHcpBasalFault, addHcpPrismFault, addHcpPy1Fault, addHcpPy1Twin, addHcpPy1TB, &
          patternFile, &
          rotate, rot, &
          verbosity, debug
  
      crystal=''
      alat=1.d0
      coa=sqrt(8./3.)   ! c/a ratio for hcp
      addHcpBasalFault=.FALSE.
      addHcpPrismFault=.FALSE.
      addHcpPy1Fault=.FALSE.
      addHcpPy1Twin=.FALSE.
      addHcpPy1TB=.FALSE.
      patternFile=''
      rotate=.false.
      rot(1:3,1:3)=ReShape( (/ 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /), (/ 3,3 /) )
      verbosity=4
      debug=.FALSE.
  
      !===========================================
      ! Read input data
      ! If input is read from "keyboard", write input in a temporary file
      IF (input_file.EQ.'-') THEN
              input_file="input.babel.temp"
              CALL WriteInput(input_file)
      END IF
              
      ! Read input data
      OPEN(file=input_file,unit=50,action='read',status='old')
      READ(50,nml=input)
      CLOSE(50)

      IF (crystal.NE.'hcp'.AND.crystal.NE.'bcc') THEN
              WRITE(0,'(2a)') 'Unknown crystal type: crystal = ', Trim(crystal)
              STOP '< PatternInit >'
      END IF

      IF ( addHcpBasalFault.AND.( crystal.NE.'hcp') ) THEN
              WRITE(0,'(a)') 'Crystal should be defined to "hcp" to add patterns&
                & corresponding to basal stacking faults'
              STOP '< PatternInit >'
      END IF
  
      IF ( addHcpPrismFault.AND.( crystal.NE.'hcp') ) THEN
              WRITE(0,'(a)') 'Crystal should be defined to "hcp" to add patterns&
                & corresponding to prismatic stacking faults'
              STOP '< PatternInit >'
      END IF
  
      IF ( addHcpPy1Fault.AND.( crystal.NE.'hcp') ) THEN
              WRITE(0,'(a)') 'Crystal should be defined to "hcp" to add patterns&
                & corresponding to pyramidal stacking faults'
              STOP '< PatternInit >'
      END IF
  
      IF ( addHcpPy1Twin.AND.( crystal.NE.'hcp') ) THEN
              WRITE(0,'(a)') 'Crystal should be defined to "hcp" to add patterns&
                & corresponding to 1st order pyramidal twins'
              STOP '< PatternInit >'
      END IF

      IF ( addHcpPy1TB.AND.( crystal.NE.'hcp') ) THEN
              WRITE(0,'(a)') 'Crystal should be defined to "hcp" to add patterns&
                & corresponding to 1st order pyramidal twin boundaries'
              STOP '< PatternInit >'
      END IF
  
      IF (verbosity.GE.verbosity_max) THEN
              WRITE(out,*)
              SELECT CASE(Trim(crystal))
              CASE('hcp')
                      WRITE(out,'(a)') 'Pattern initialization for hcp crystal'
                      WRITE(out,'(a,g14.6)') '  lattice parameter :  alat = ', alat
                      WRITE(out,'(a,g14.6)') '  c/a ratio :          coa  = ', coa
                      IF (addHcpBasalFault) THEN
                              WRITE(out,'(a)') 'patterns corresponding to basal stacking faults will be added'
                      END IF
                      IF (addHcpPrismFault) THEN
                              WRITE(out,'(a)') 'patterns corresponding to prismatic stacking faults will be added'
                              WRITE(out,'(a)') '  fault plane: (10-10)'
                              WRITE(out,'(a)') '  fault vector: 1/6 [1-210]'
                      END IF
                      IF (addHcpPy1Fault) THEN
                              WRITE(out,'(a)') 'patterns corresponding to pyramidal 1 stacking faults will be added'
                              WRITE(out,'(a)') '  fault planes: (10-11) and (-1011)'
                              WRITE(out,'(a)') '  fault vector: 1/6 [1-210] (plus a small orthogonal component)'
                      END IF
                      IF (addHcpPy1Twin) THEN
                              WRITE(out,'(a)') 'patterns corresponding to 1st order pyramidal twins will be added'
                      END IF
                      IF (addHcpPy1TB) THEN
                              WRITE(out,'(a)') 'patterns corresponding to 1st order pyramidal twin boundaries will be added'
                      END IF
              CASE('bcc')
                      WRITE(out,'(a)') 'Pattern initialization for bcc crystal'
                      WRITE(out,'(a,g14.6)') '  lattice parameter :  alat = ', alat
              CASE DEFAULT
                      WRITE(out,'(2a)') 'Unknown crystal type: crystal = ', Trim(crystal)
              END SELECT
              IF (rotate) THEN
                      WRITE(out,'(a)') ' Rotation of the patterns according to the matrix'
                      WRITE(out,'(a,3g14.6,a)') '   | ', rot(1,1:3), ' |'
                      WRITE(out,'(a,3g14.6,a)') '   | ', rot(2,1:3), ' |'
                      WRITE(out,'(a,3g14.6,a)') '   | ', rot(3,1:3), ' |'
                      WRITE(out,*)
              END IF
      END IF
    END SUBROUTINE read_inputPattern

END PROGRAM PatternInit
