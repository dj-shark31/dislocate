PROGRAM PatternDetect

  USE babel_data
  USE structure_module
  USE neighboursModule
  USE gradDisplacementModule
  USE gradElasticDisplacementModule
  USE strainFromDisplacementModule
  USE nyeTensorModule 

  IMPLICIT NONE
 
  INTEGER, parameter :: verbosity_max=1

  ! Number of neighbours for each atom (reference structure)
  INTEGER, dimension(:), allocatable :: nNeigh0
  ! Neighbour indexes for each atom (reference structure)
  INTEGER, dimension(:,:), allocatable :: iNeigh0

  INTEGER :: nAux_real, nAux_int, n, i
  CHARACTER(len=50), dimension(:), allocatable :: aux_title
  REAL(kind(0.d0)), dimension(:,:), allocatable :: aux_real
  INTEGER, dimension(:,:), allocatable :: aux_int  
  !INTEGER, EXTERNAL :: iArgc
  CHARACTER(len=9) :: fileType

  ! Name of the input file defining the dislocation
  CHARACTER(len=100) :: input_file

  program_name="patternDetect"

  ! Read name of the input file
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF

  ! Read input parameters in file 'input.dat'
  CALL read_inputPattern(6)

  ! Rotation
  IF (rotate) THEN  
          ! Rotate atom cartesian coordinates
          IF (verbosity.GE.verbosity_max) THEN
                  WRITE(6,'(a)') ' Rotation applied'
                  WRITE(6,'(a,3g14.6,a)') '   | ', rot(1,1:3), ' |'
                  WRITE(6,'(a,3g14.6,a)') '   | ', rot(2,1:3), ' |'
                  WRITE(6,'(a,3g14.6,a)') '   | ', rot(3,1:3), ' |'
                  WRITE(6,*)
          END IF
          CALL rotate_cell(rot, xp0, im0, at0, 6, verbose=verbosity.GE.3) 
          rotate=.FALSE.
  END IF

  ! Duplicate unit cell
  IF (duplicate) THEN
          CALL duplicate_cell(lat, xp0, iTyp0, im0, at0, 6, verbose=.FALSE.)
          duplicate=.FALSE.
          Lat(1:3)=1
  END IF

  ! Add solid displacement
  IF (translate) THEN
          DO i=1, im
            xp(1:3,i) = xp(1:3,i) + uTranslate(1:3)
            xp0(1:3,i) = xp0(1:3,i) + uTranslate(1:3)
          END DO
  END IF

  ! Clip atoms in unit cell applying periodic boundary conditions
  IF (clipAtom) THEN
          CALL Clip_Atoms(xp0, im0, at0, 6, verbosity.GE.3)
          clipAtom=.false.
  END IF

  ! Calculate neighbours table for reference structure
  CALL InitNeighbours(imm)
  ALLOCATE(nNeigh0(1:imm))               ; nNeigh0(:) = 0
  ALLOCATE(iNeigh0(1:max_nNeigh, 1:imm)) ; iNeigh0(:,:) = 0
  CALL BuildNeighbours(im0, xp0, at0, nNeigh0, iNeigh0, 6)

  ! Search different patterns in reference structure file
  CALL InitRefPattern(im0, xp0, at0, nNeigh0, iNeigh0, 6) 

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

  ! Properties printed on output
  nAux_real=0 ; nAux_int=0; n=0
  IF (out_neighbours)              nAux_int  = nAux_int  + 1
  IF (out_pattern)                 nAux_int  = nAux_int  + 1

  IF ( (nAux_int.NE.0) .OR. (nAux_real.NE.0) ) THEN

          ALLOCATE(aux_title(1:nAux_real+nAux_int))
          IF (nAux_int.NE.0) THEN
                  ALLOCATE(aux_int(1:nAux_int,1:imm))
                  aux_int(:,:)=0
          END IF
          IF (nAux_real.NE.0) THEN
                  ALLOCATE(aux_real(1:nAux_real,1:imm)) 
                  aux_real(:,:) = 0.d0
          END IF

          IF (out_neighbours) THEN      ! Print number of neighbours for each atom
                  aux_int(n+1,:) = nNeigh0(:)
                  aux_title(n+1)='Number of neighbours'
                  n=n+1
          END IF
          IF (out_pattern) THEN      ! Print number of neighbours for each atom
                  aux_int(n+1,:) = refPattern(:)
                  aux_title(n+1)='Pattern index'
                  n=n+1
          END IF
  END IF

  IF (n.NE.nAux_real+nAux_int) THEN
          WRITE(0,'(a,i0)') ' n         = ', n
          WRITE(0,'(a,i0)') ' nAux_int = ', nAux_int
          WRITE(0,'(a,i0)') ' nAux_real = ', nAux_real
  END IF

  ! ==== Output structure ========================
  IF (verbosity.GE.verbosity_max) THEN
          WRITE(6,*)
          IF (outXyz)       WRITE(6,'(2a)') 'Write output structure with Xyz format in file ', outFile
          IF (outOnlyAtoms) WRITE(6,'(2a)') 'Write one line per atom in file ', outFile
          IF (outCfg)       WRITE(6,'(2a)') 'Write output structure with Cfg format in file ', outFile
          IF (outLammps)    WRITE(6,'(2a)') 'Write output structure with Lammps dump format in file ', outFile
          IF (outXyz) THEN
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': atom label'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': x coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': y coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 4, ': z coordinate (A)'
                  DO n=1, nAux_int+nAux_real
                     WRITE(6,'(a,i2,2a)') '   column ', n+4, ': ', aux_title(n)
                  END DO
          END IF
          IF (outOnlyAtoms) THEN
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': x coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': y coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': z coordinate (A)'
                  DO n=1, nAux_int+nAux_real
                     WRITE(6,'(a,i2,2a)') '   column ', n+3, ': ', aux_title(n)
                  END DO
          END IF
          IF (outCfg) THEN
                  DO n=1, nAux_int+nAux_real
                     WRITE(6,'(a,i2,2a)') '   property ', n-1, ': ', aux_title(n)
                  END DO
          END IF
          IF (outLammps) THEN
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': atom id'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': atom type'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': reduced coordinate xsu'
                  WRITE(6,'(a,i2,a)') '   column ', 4, ':                    ysu'
                  WRITE(6,'(a,i2,a)') '   column ', 5, ':                    zsu'
                  DO n=1, nAux_int+nAux_real
                     WRITE(6,'(a,i2,2a)') '   column ', n+5, ': ', aux_title(n)
                  END DO
          END IF
          WRITE(6,*)
  END IF

  IF ( outXyz.OR.outOnlyAtoms.OR.outCfg.OR.outGin.OR.outLisa.OR.outSiesta.OR.outNDM.OR.outLammps ) THEN

          IF (outXyz)       fileType='xyz'
          IF (outOnlyAtoms) fileType='onlyAtoms'
          IF (outCfg)       fileType='cfg'
          IF (outLammps)    fileType='lammps'
          IF (outGin)    WRITE(0,'(a)') "Gin type not supported for output structure file"
          IF (outLisa)   WRITE(0,'(a)') "Lisa type not supported for output structure file"
          IF (outSiesta) WRITE(0,'(a)') "Siesta type not supported for output structure file"
          IF (outNDM)    WRITE(0,'(a)') "Binary NDM type not supported for output structure file"
          IF ( (nAux_int.GT.0).AND.(nAux_real.GT.0) ) THEN
                  CALL WriteStructure(xp0, iTyp0, im0, at0, outFile, fileType, &
                        nAux_int=nAux_int, aux_int=aux_int(:,:), &
                        nAux_real=nAux_real, aux_real=aux_real(:,:), &
                        aux_title=aux_title(:))
          ELSE IF (nAux_int.GT.0) THEN
                  CALL WriteStructure(xp0, iTyp0, im0, at0, outFile, fileType, &
                        nAux_int=nAux_int, aux_int=aux_int(:,:), &
                        aux_title=aux_title(:))
          ELSE IF (nAux_real.GT.0) THEN
                  CALL WriteStructure(xp0, iTyp0, im0, at0, outFile, fileType, &
                        nAux_real=nAux_real, aux_real=aux_real(:,:), &
                        aux_title=aux_title(:))
          ELSE
                  CALL WriteStructure(xp0, iTyp0, im0, at0, outFile, fileType)
          END IF

  END IF

CONTAINS

  SUBROUTINE read_inputPattern(out)
  
      USE babel_data
      USE structure_module
      USE neighboursModule, ONLY : rNeigh, max_nNeigh
      USE GradElasticDisplacementModule, ONLY : patternFile
      IMPLICIT NONE
      INTEGER, intent(in) :: out
  
      REAL(kind(0.d0)) :: alat     ! Lattice parameter (A)
  
      NAMELIST /input/ alat, &
          translate, uTranslate, &
          clipAtom, &
          imm, duplicate, lat, rotate, rot, &
          refXyz, refCfg, refGin, refLisa, refSiesta, refNDM, refFile, refLammps, refPoscar, &
          patternFile, &
          verbosity, debug, &
          nTypes, mass, label, &
          rNeigh, max_nNeigh, out_neighbours, &
          out_pattern, &
          out_alat, outXyz, outCfg, outOnlyAtoms, outLammps, outFile
  
      CHARACTER(len=6) :: fileType
  
      !===========================================
      ! Initialisation
      CALL Random_Seed()
  
      refFile='-'
      refXyz=.FALSE. ; refCfg=.FALSE. ; refGin=.FALSE. ; refLisa=.FALSE. ; refSiesta=.FALSE. ; refNDM=.FALSE. ; refLammps=.FALSE. ; refPoscar=.FALSE.
      imm=0 
      duplicate=.FALSE.
      lat(1:3)=1
      rotate=.FALSE.
      rot(1:3,1:3)=matId(1:3,1:3) ! Identity matrix
      translate=.FALSE.
      uTranslate(1:3) = 0.d0
      alat=1.d0
      patternFile=''
      rNeigh=-1.d0
      max_nNeigh=16
      outFile='-' 
      outXyz=.FALSE. ; outCfg=.FALSE. ; outGin=.FALSE. ; outLisa=.FALSE. ; outSiesta=.FALSE. ; outNDM=.FALSE. ; outOnlyAtoms=.FALSE. ; outLammps=.FALSE.
      out_neighbours=.FALSE.
      out_pattern=.FALSE.
      clipAtom=.FALSE.
      verbosity=4
      debug=.FALSE.
  
      CALL InitStructure()        ! Initialize out_alat, nTypes, mass(:), label(:)
  
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
  
      IF (rNeigh.GT.0.d0) rNeigh=rNeigh*alat
  
      !===========================================
      ! Read reference structure file
      IF ( Count( (/refXyz,refCfg,refGin,refLisa,refSiesta,refNDM,refLammps,refPoscar/) ).GT.1 ) THEN
              WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin,  Lisa, Siesta, Ndm, Lammps, Poscar formats for reference structure'
              STOP '< PatternDetect >'
      END IF
      IF (refXyz.OR.refCfg.OR.refGin.OR.refLisa.OR.refSiesta.OR.refNDM.OR.refLammps.OR.refPoscar) THEN
              ! Reference file is an input file for program displacement
              !   for babel, it is an output file
              IF (refXyz)    fileType='xyz'
              IF (refCfg)    fileType='cfg'
              IF (refGin)    fileType='gin'
              IF (refSiesta) fileType='siesta'
              IF (refLisa)   fileType='lisa'
              IF (refNDM)    fileType='ndm'
              IF (refLammps) fileType='lammps'
              IF (refPoscar) fileType='poscar'
              CALL ReadStructure(xp0, iTyp0, im0, imm, at0, refFile, fileType)
              IF (Allocated(keep)) Deallocate(keep)
              ALLOCATE(keep(1:imm))
              keep(:)=.TRUE.
      END IF
  
  
      !===========================================
      IF (translate) THEN
              uTranslate(:) = alat*uTranslate(:)
              IF (Sum(uTranslate(:)**2).LE.1.d-12) translate=.FALSE.
      END IF
  
      !===========================================
      IF ( Count( (/outXyz,outCfg,outOnlyAtoms,outLammps/) ).GT.1 ) THEN
              WRITE(0,'(a)') 'Choose between Xyz, OnlyAtoms and Cfg formats for output structure'
              STOP '< Read_Babel >'
      END IF
  
    END SUBROUTINE read_inputPattern

END PROGRAM PatternDetect
