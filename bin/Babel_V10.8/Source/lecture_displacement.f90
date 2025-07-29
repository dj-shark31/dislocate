MODULE DisplacementModule

  SAVE

  ! Name of the input file defining the dislocation
  CHARACTER(len=100) :: input_file

  INTEGER, parameter, private :: verbosity_max=2

CONTAINS

SUBROUTINE Read_Displacement(out)

    USE babel_data
    USE structure_module
    USE neighboursModule, ONLY : rNeigh, max_nNeigh
    USE GradElasticDisplacementModule, ONLY : patternAngleThreshold, &
            patternFile, patternSelectionMethod
    USE StrainFromDisplacementModule, ONLY : EulerLagrange
    IMPLICIT NONE
    INTEGER, intent(in) :: out

    REAL(kind(0.d0)) :: alat     ! Lattice parameter (A)
    REAL(kind(0.d0)) :: at_norm2, inv_norm
    REAL(kind(0.d0)), dimension(1:3,1:3) :: at_temp
    REAL(kind(0.d0)), dimension(1:3) :: ex, ey, ez
    INTEGER :: i, n
    CHARACTER(len=6) :: fileType

    ! Format of the input structure files
    LOGICAL :: inpXyz, inpCfg, inpGin, inpLisa, inpSiesta, inpNDM, inpLammps, inpPoscar

    NAMELIST /input/ alat, &
        translate, uTranslate, fixGravity, &
        clipDisplacement, clipAtom, &
        imm, duplicate, lat, rotate, rot, at, &
        inpXyz, inpCfg, inpGin, inpLisa, inpSiesta, inpNDM, inpLammps, inpPoscar, inpFile, &
        refXyz, refCfg, refGin, refLisa, refSiesta, refNDM, refLammps, refPoscar, refFile, &
        verbosity, debug, &
        nTypes, mass, label, &
        initial, &
        rNeigh, max_nNeigh, out_neighbours, &
        patternFile, patternSelectionMethod, &
        patternAngleThreshold, out_pattern, &
        out_nye, out_BurgersDensity, bNye, lNye, &
        out_strain, out_elasticStrain, out_plasticStrain, EulerLagrange, &
        out_gradDisplacement, out_gradElasticDisplacement, out_displacement, &
        out_alat, outXyz, outCfg, outOnlyAtoms, outLammps, outFile

    !===========================================
    ! Initialisation
    CALL Random_Seed()

    inpFile='-' ; refFile='-'
    inpXyz=.FALSE. ; inpCfg=.FALSE. ; inpGin=.FALSE. ; inpLisa=.FALSE. ; inpSiesta=.FALSE. ; inpNDM=.FALSE. ; inpLammps=.FALSE. ; inpPoscar=.FALSE.
    refXyz=.FALSE. ; refCfg=.FALSE. ; refGin=.FALSE. ; refLisa=.FALSE. ; refSiesta=.FALSE. ; refNDM=.FALSE. ; refLammps=.FALSE. ; refPoscar=.FALSE.
    imm=0
    duplicate=.FALSE.
    lat(1:3)=1
    rotate=.FALSE.
    rot(1:3,1:3)=matId(1:3,1:3) ! Identity matrix
    translate=.FALSE.
    uTranslate(1:3) = 0.d0
    clipAtom=.FALSE.
    clipDisplacement=.FALSE.
    alat=1.d0
    fixGravity=.FALSE.
    patternFile=''
    patternAngleThreshold=10.d0
    patternSelectionMethod=3
    rNeigh=-1.d0
    max_nNeigh=16
    bNye(:)=0.d0
    lNye=0.d0
    outFile='-'
    outXyz=.FALSE. ; outCfg=.FALSE. ; outGin=.FALSE. ; outLisa=.FALSE. ; outSiesta=.FALSE. ; outNDM=.FALSE. ; outOnlyAtoms=.FALSE. ; outLammps=.FALSE.
    initial=.FALSE.
    out_displacement=.true.
    out_gradDisplacement=.FALSE.
    out_strain=.FALSE.
    EulerLagrange=.FALSE.
    out_gradElasticDisplacement=.FALSE.
    out_elasticStrain=.FALSE.
    out_plasticStrain=.FALSE.
    out_nye=.FALSE.
    out_BurgersDensity=.FALSE.
    out_neighbours=.FALSE.
    out_pattern=.FALSE.
    at(1:3,1:3)=0.d0
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
    ! Basis vectors for input structure could be defined here or read in structure file
    at_norm2=SUM( at(1:3,1:3)**2 )
    IF (at_norm2.GT.Distance_Zero2) THEN
            ! Basis vectors defined
            at(:,:) = alat*at(:,:)
            at_defined=.TRUE.
    ELSE
            at_defined=.FALSE.
    END IF


    !===========================================
    ! Read input structure file

    IF ( Count( (/inpXyz,inpCfg,inpGin,inpLisa,inpSiesta,inpNDM,inpLammps,inpPoscar/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin, Lisa,Siesta, Lammps, and PosCar formats for input structure'
            STOP '< Read_Displacement >'
    END IF
    IF (inpXyz.OR.inpCfg.OR.inpGin.OR.inpLisa.OR.inpSiesta.OR.inpNDM.OR.inpLammps.OR.inpPoscar) THEN
            IF (inpXyz)    fileType='xyz'
            IF (inpCfg)    fileType='cfg'
            IF (inpGin)    fileType='gin'
            IF (inpSiesta) fileType='siesta'
            IF (inpLisa)   fileType='lisa'
            IF (inpNDM)    fileType='ndm'
            IF (inpLammps) fileType='lammps'
            IF (inpPoscar) fileType='poscar'
            CALL ReadStructure(xp, iTyp, im, imm, at_temp, inpFile, fileType)
            IF (Allocated(keep)) Deallocate(keep)
            ALLOCATE(keep(1:imm))
            keep(:)=.TRUE.
    END IF

    !===========================================
    ! Basis vectors for input structure could be defined here or read in structure file
    at_norm2=SUM( at_temp(1:3,1:3)**2 )
    IF (at_norm2.GT.Distance_Zero2) THEN
            IF (at_defined) THEN
                    WRITE(0,'(a)') 'Periodicity vectors at(1:3,i)&
                        & defined in input and structure files '
                    STOP '< Read_Displacement >'
            ELSE
                    at(:,:) = at_temp
                    at_defined=.TRUE.
            END IF
    END IF

    !===========================================
    ! Basis vectors for reference structure could be defined here or read in structure file
    at_norm2=SUM( at0(1:3,1:3)**2 )
    IF (at_norm2.GT.Distance_Zero2) THEN
            ! Basis vectors defined
            at0(:,:) = alat*at0(:,:)
            at0_defined=.TRUE.
    ELSE
            at0_defined=.FALSE.
    END IF

    !===========================================
    ! Read reference structure file

    IF ( Count( (/refXyz,refCfg,refGin,refLisa,refSiesta,refNDM,refLammps,refPoscar/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin,  Lisa, Siesta, Lammps, and Poscar formats for reference structure'
            STOP '< Read_Displacement >'
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
            CALL ReadStructure(xp0, iTyp0, im0, imm, at_temp, refFile, fileType)
            IF (Allocated(keep)) Deallocate(keep)
            ALLOCATE(keep(1:imm))
            keep(:)=.TRUE.

            !===========================================
            ! Basis vectors for input structure could be defined here or read in structure file
            at_norm2=SUM( at_temp(1:3,1:3)**2 )
            IF (at_norm2.GT.Distance_Zero2) THEN
                    IF (at0_defined) THEN
                            WRITE(0,'(a)') 'Periodicity vectors at0(1:3,i)&
                                & defined in input and structure files '
                            STOP '< Read_Displacement >'
                    ELSE
                            at0(:,:) = at_temp
                            at0_defined=.TRUE.
                    END IF
            ELSE
                    at0 = at
                    at0_defined = at_defined
            END IF

            ! Check number of atoms
            IF (im0.NE.im) THEN
                    WRITE(0,'(a)') 'Both input and reference structure files do not contain the&
                          & same number of atoms'
                    WRITE(0,'(a,i0)') '  im0 = ', im0
                    WRITE(0,'(a,i0)') '  im  = ', im 
                    STOP '< Read_Displacement >'
            END IF

           ! Check that all atoms are of the same type
           DO i=1, im
              IF (iTyp0(i).nE.iTyp(i)) THEN
                      WRITE(0,'(a,i0,a)') 'Atom ', i, ' is not of the same type in both input and reference structure files'
                      WRITE(0,'(a,i0)') '  type in input file: ', iTyp(i)
                      WRITE(0,'(a,i0)') '  type in reference file: ', iTyp0(i)
                      STOP '< Read_Displacement >'
              END IF
           END DO
    END IF


    !===========================================
    IF (translate) THEN
            uTranslate(:) = alat*uTranslate(:)
            IF (Sum(uTranslate(:)**2).LE.1.d-12) translate=.FALSE.
    END IF

    !===========================================
    IF ( Count( (/outXyz,outCfg,outOnlyAtoms,outLammps/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, OnlyAtoms and Cfg formats for output structure'
            STOP '< Read_Displacement >'
    END IF

END SUBROUTINE read_displacement

END MODULE displacementModule
