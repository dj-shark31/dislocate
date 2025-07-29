SUBROUTINE Lecture_input(inp, alat)

    USE babel_data
    USE neighboursModule, ONLY : rNeigh, max_nNeigh
    USE LineCouple_elasticity_ani
    USE LineForce_elasticity_ani
    USE structure_module
    USE xyz_module
    USE math
    USE elasticity_ani
    USE slabModule
    USE symmetryXmlModule
    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    REAL(kind(0.d0)), intent(inout) :: alat     ! Lattice parameter (A)

    ! Elastic constants (cubic symmetry)
    LOGICAL :: cubic_elasticity, hexagonal_elasticity
    REAL(kind(0.d0)) :: C11, C12, C13, C33, C44
    REAL(kind(0.d0)), dimension(1:3,1:3) :: sStrain
    REAL(kind(0.d0)), dimension(1:6) :: sVoigt, eVoigt
    REAL(kind(0.d0)) :: epsi, tau, volume
    CHARACTER(len=100) :: displacementFile
    REAL(kind(0.d0)) :: displacementFileFactor
    REAL(kind(0.d0)) :: ex_norm2, inv_norm
    REAL(kind(0.d0)), dimension(1:3,1:3) :: at_temp
    REAL(kind(0.d0)), dimension(1:3) :: ex, ey, ez
    INTEGER :: i, n
    CHARACTER(len=6) :: fileType
    TYPE(sym3D_group_t) :: symGroup0

    ! Format of the input structure files
    LOGICAL :: inpXyz, inpCfg, inpGin, inpLisa, inpSiesta, inpNDM, inpLammps, inpPoscar
    
    ! Xml file containing definition of symmetry operations
    CHARACTER(len=100) :: symFile

    ! Temporary variables needed to read displacement file
    INTEGER :: nTypes_temp, im_temp
    CHARACTER(len=5), dimension(1:max_nTypes) :: label_temp
    REAL(kind(0.d0)), dimension(:,:), allocatable :: xp_temp
    INTEGER, dimension(:), allocatable :: iTyp_temp
    LOGICAL :: add_cut


    NAMELIST /input/  &
        strain, epsi, eStrain, tau, sStrain, induced_homogeneous_strain, &
        translate, uTranslate, fixGravity, &
        imm, duplicate, lat, rotate, rot, symmetrize, symFile, &
        alat, &
        C11, C12, C13, C33, C44, CVoigt, CVoigt_noise, &
        cubic_elasticity, hexagonal_elasticity, anisotropic_elasticity, &
        xImages, yImages, zImages, nxImages, nyImages, nzImages, &
        clipDisplacement, clipAtom, remove_cut, add_cut, &
        max_Euler, delta_Euler, &
        xNoise, &
        rc, factorE, &
        inpXyz, inpCfg, inpGin, inpLisa, inpSiesta, inpNDM, inpLammps, inpPoscar, inpFile, &
        nTypes, mass, label, at, &
        outXyz, outOnlyAtoms, outCfg, outGin, outLisa, outSiesta, outLammps, outPoscar, outFile, & 
        out_alat, initial, &
        out_neighbours, rNeigh, max_nNeigh, &
        out_id, &
        out_displacement, out_elasticStrain, out_stress, out_pressure, out_VonMises, &
        out_core, &
        out_Ebinding, pImpurity, &
        displacementFile, displacementFileFactor, &
        !&
        refXyz, refCfg, refGin, refLisa, refSiesta, refNDM, refLammps, refPoscar, refFile, &  ! OBSOLETE (or soon)
        out_x, out_y, out_z, &
        verbosity, debug, &
        !&
        slab, nxSlab, nySlab, nzSlab, & ! OBSOLETE (or not working)
        LineForce, lLineForce, fLineForce, cLineForce, cutLineForce, nlf ! OBSOLETE (or soon)

    !===========================================
    ! Initialisation
    CALL Random_Seed()

    strain=.FALSE.
    epsi = 1.d0
    eStrain(1:3,1:3) = 0.d0
    tau = 1.d0
    sStrain(1:3,1:3) = 0.d0
    induced_homogeneous_strain=.true.
    translate=.FALSE.
    uTranslate(1:3) = 0.d0
    fixGravity=.FALSE.
    imm=0
    duplicate=.FALSE.
    lat(1:3)=1
    rotate=.FALSE.
    rot(1:3,1:3)=matId(1:3,1:3) ! Identity matrix
    symmetrize=.FALSE.
    symFile=''
    C11 = 0.d0 ; C12 = 0.d0 ; C13 = 0.d0 ; C33 = 0.d0 ; C44 = 0.d0
    CVoigt(:,:)=0.d0
    CVoigt_noise=0.d0
    inv_CVoigt(:,:)=0.d0
    cubic_elasticity=.FALSE. 
    hexagonal_elasticity=.FALSE.
    anisotropic_elasticity=.FALSE. 
    xImages=.FALSE.
    yImages=.FALSE.
    zImages=.FALSE.
    nxImages=10
    nyImages=10
    nzImages=10
    clipDisplacement=.FALSE.
    clipAtom=.FALSE.
    remove_cut=.false.
    add_cut=.false.
    max_Euler=0 ! default = Lagrangi, atan coordinates
    delta_Euler=1.d-5
    xNoise=0.d0        ! Noise to add to atomic positions
    rc=1.d0
    factorE=1.d9*1.d-30/1.602176530000000045d-19        ! GPa.A^3 => eV
    inpFile='-' ; outFile='-' ; refFile='-'
    inpXyz=.FALSE. ; inpCfg=.FALSE. ; inpGin=.FALSE. ; inpLisa=.FALSE. ; inpSiesta=.FALSE. ; inpNDM=.FALSE. ; inpLammps=.FALSE. ; inpPoscar=.FALSE.
    outXyz=.FALSE. ; outCfg=.FALSE. ; outGin=.FALSE. ; outLisa=.FALSE. ; outSiesta=.FALSE. ; outNDM=.FALSE. ; outOnlyAtoms=.FALSE. ; outLammps=.FALSE. ; outPoscar=.FALSE.
    refXyz=.FALSE. ; refCfg=.FALSE. ; refGin=.FALSE. ; refLisa=.FALSE. ; refSiesta=.FALSE. ; refNDM=.FALSE. ; refLammps=.FALSE. ; refPoscar=.FALSE.
    at(1:3,1:3)=0.d0
    initial=.FALSE.
    out_neighbours=.FALSE.
    rNeigh=-1.d0
    max_nNeigh=16
    out_id=.FALSE.
    out_displacement=.FALSE.
    out_elasticStrain=.FALSE.
    out_stress=.FALSE.
    out_pressure=.FALSE.
    out_VonMises=.FALSE.
    out_core=.false.
    out_Ebinding=.FALSE.
    pImpurity(1:3,1:3) = 0.d0
    out_x=.FALSE.
    out_y=.FALSE.
    out_z=.FALSE.
    verbosity=4
    debug=.FALSE.
    read_uFile=.FALSE.
    displacementFile=''
    displacementFileFactor=1.

    LineForce=.FALSE.
    nlf=0
    lLineForce(1:3,:) = 0.d0
    fLineForce(1:3,:) = 0.d0
    cLineForce(1:3,:) = 0.d0
    cutLineForce(1:3,:) = 0.d0

    slab=.false.
    nxSlab=1
    nySlab=1
    nzSlab=1
    CALL InitStructure()        ! Initialize out_alat, nTypes, mass(:), label(:)


    !===========================================
    ! Read input data
    READ(inp,nml=input)
    remove_cut=remove_cut.OR.add_cut

    IF (program_name.EQ."babel") THEN
            IF (clipDisplacement) THEN
                    WRITE(0,'(a)') 'Option clipDisplacement does not work with program Babel'
                    WRITE(0,'(a)') 'You need to remove it from your input file'
                    STOP '< Lecture_input >'
            END IF
    END IF

    rc=rc*alat
    IF (rNeigh.GT.0.d0) rNeigh=rNeigh*alat
    IF (xNoise.GT.0.d0) xNoise=xNoise*alat

    !===========================================
    ! Basis vectors could be defined here or read in structure file
    IF (SUM( at(1:3,1:3)**2 ).GT.Distance_Zero2) THEN
            ! Basis vectors defined
            at(:,:) = alat*at(:,:)
            at_defined=.TRUE.
    ELSE
            at_defined=.FALSE.
    END IF


    !===========================================
    ! Read input structure file

    IF ( Count( (/inpXyz,inpCfg,inpGin,inpLisa,inpSiesta,inpNDM,inpLammps,inpPoscar/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin, Lisa, Siesta, Lammps, and Poscar formats for input structure'
            STOP '< Lecture_input >'
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
    ! Basis vectors could be defined here or read in structure file
    IF (SUM( at_temp(1:3,1:3)**2 ).GT.Distance_Zero2) THEN
            IF (at_defined) THEN
                    WRITE(0,'(a)') 'Periodicity vectors at(1:3,i)&
                        & defined in input and structure files '
                    STOP '< Lecture_input >'
            ELSE
                    at(:,:) = at_temp
                    at_defined=.TRUE.
            END IF
    END IF

    ! Check that basis vectors are linearly independant
    IF (at_defined) THEN
            volume=abs(MatDet(at))
            IF (abs(volume).LE.distance_Zero**3) THEN
                    WRITE(0,'(a)') 'Periodicity vectors are not linearly indepedant'
                    WRITE(0,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
                    WRITE(0,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
                    WRITE(0,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
                    WRITE(0,'(a,g20.12)')  '       corresponding volume: ', volume
                    STOP '< Lecture_input >'
            END IF
    END IF

    !===========================================
    ! Read reference structure file

    IF ( Count( (/refXyz,refCfg,refGin,refLisa,refSiesta,refNDM,refLammps,refPoscar/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin,  Lisa, and Siesta formats for reference structure'
            STOP '< Lecture_input >'
    END IF
    IF ( (program_name.EQ."displacement") .AND. &
        (refXyz.OR.refCfg.OR.refGin.OR.refLisa.OR.refSiesta) ) THEN
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

            ! Periodicity vectors
            IF (SUM( at0(1:3,1:3)**2 ).LE.Distance_Zero2) at0=at

            ! Check number of atoms
            IF (im0.NE.im) THEN
                    WRITE(0,'(a)') 'Both input and reference structure files do not contain the&
                          & same number of atoms'
                    WRITE(0,'(a,i0)') '  im0 = ', im0
                    WRITE(0,'(a,i0)') '  im  = ', im 
                    STOP '< Lecture_input >'
            END IF

           ! Check that all atoms are of the same type
           DO i=1, im
              IF (iTyp0(i).nE.iTyp(i)) THEN
                      WRITE(0,'(a,i0,a)') 'Atom ', i, ' is not of the same type in both input and reference structure files'
                      WRITE(0,'(a,i0)') '  type in input file: ', iTyp(i)
                      WRITE(0,'(a,i0)') '  type in reference file: ', iTyp0(i)
                      STOP '< Read_Babel >'
              END IF
           END DO
    END IF

    !===========================================
    ! Read displacement file

    IF (displacementFile.NE.'') THEN
            read_uFile=.TRUE.
            IF (Allocated(uFile)) Deallocate(uFile)
            ALLOCATE(uFile(1:3,1:imm))
            uFile(:,:)=0.d0
            OPEN(file=displacementFile,unit=50,status='old',action='read')
            ALLOCATE(xp_temp(1:3,1:imm))
            ALLOCATE(iTyp_temp(1:imm))
            nTypes_temp=0
            CALL ReadXyz(xp_temp, iTyp_temp, im_temp, at_temp, nTypes_temp, label_temp, 50, nAux=3, aux=uFile)
            DEALLOCATE(xp_temp)
            DEALLOCATE(iTyp_temp)
            CLOSE(50)
            uFile(:,:) = displacementFileFactor*uFile(:,:)
    ELSE
            read_uFile=.FALSE.
    END IF

    !===========================================
    IF (symmetrize) THEN
            IF ( symFile.EQ.'' ) THEN
                    WRITE(0,'(a)') 'symmetrize=.true.'
                    WRITE(0,'(a)') 'You need to define variable symFile with the xml file containing symmetry operations'
                    STOP '< Lecture_input >'
            ELSE
                    CALL ReadSymmetryFileXml(symFile, symGroup0)
                    CALL BuildSpaceGroup3D(symGroup0, symGroup)
                    IF (at_defined) THEN 
                            CALL CheckBasisSym3Dgroup(at, symGroup, 6)
                    END IF
            END IF
    END IF
    !===========================================
    ! Periodic boundary conditions
    ! number of dislocation images to consider in each direction
    IF (.NOT.xImages) nxImages=0
    IF (.NOT.yImages) nyImages=0
    IF (.NOT.zImages) nzImages=0

    IF ( (xImages.OR.yImages.OR.zImages).AND.(.NOT.at_defined) ) THEN
            WRITE(0,'(a)') 'You need to define periodicity vectors at(1:3,i) &
                &to use periodic boundary conditions'
            STOP '< Lecture_input >'
    END IF

    !===========================================
    ! Initialization for line-force 
    IF (LineForce) THEN

            IF (nlf.LE.0) THEN
                    WRITE(0,'(a)') 'You need to enter a positive number of line-forces (nlf)'
                    STOP '< Lecture_input >'
            ELSEIF (nlf.GT.max_nlf) THEN
                    WRITE(0,'(a,i0)') 'Maximal number of line-forces (max_nlf): ', max_nlf
                    WRITE(0,'(a,i0)') 'Current number of line-forces (nlf): ', nlf
                    WRITE(0,'(a)') 'you need to change value in data.f90 and to re-compile'
                    STOP '< Lecture_input >'
            END IF
                    
            IF (.NOT.( cubic_elasticity.OR.hexagonal_elasticity.OR.anisotropic_elasticity)) THEN
                    WRITE(0,'(a)') 'Cubic_Elasticity, Hexagonal_Elasticity, or Anistropic_Elasticity &
                        &need to be set to .TRUE. for &
                        &line-force calculation'
                    STOP '< Lecture_input >'
            END IF

            loop_force1: DO n=1, nlf    ! Loop on all line-forces

                    ! Check cutting direction
                    IF ( Abs( Sum( cutLineForce(1:3,n)*lLineForce(1:3,n) ) ).GT.Distance_Zero ) THEN
                            WRITE(0,'(a,i0)') 'Line-force ', n
                            WRITE(0,'(a)') 'Cutting direction cutLineForce(1:3) has to be &
                                &orthogonal to the line direction'
                            STOP '< Lecture_input >'
                    END IF

                    ! Multiply center of line-force by lattice parameter
                    cLineForce(1:3,n) = alat*cLineForce(1:3,n)

                    ! Normalize line-force direction
                    inv_norm = 1.d0/Sqrt( Sum( lLineForce(1:3,n)**2 ) )
                    lLineForce(1:3,n) = inv_norm*lLineForce(1:3,n)
                    ez(1:3) = lLineForce(1:3,n)

                    ! Orientate the crystal
                    !  z: line-force direction
                    !  y: 
                    !  x: cutting direction
                    ex_norm2 = Sum( cutLineForce(1:3,n)**2 )
                    IF (ex_norm2.GT.Distance_Zero2) THEN
                            ! Cutting direction is defined
                            inv_norm = 1.d0/Sqrt( ex_norm2 )
                            ex(1:3) = inv_norm*cutLineForce(1:3,n)
                    ELSE
                            ! Try vector perpendicular to lLineForce(1:3) and (1,0,0)
                            ex(1)=0.d0 ; ex(2)=lLineForce(3,n) ; ex(3)=-lLineForce(2,n)
                            ex_norm2 = Sum( ex(1:3)**2 )
                            IF (ex_norm2.LE.Distance_Zero2) THEN
                                    ! Try vector perpendicular to lLineForce(1:3) and (0,1,0)
                                    ex(1)=-lLineForce(3,n) ; ex(2)=0.d0 ; ex(3)=lLineForce(1,n)
                            ENDIF
                            inv_norm = 1.d0/Sqrt( Sum( ex(1:3)**2 ) )
                            ex(1:3) = inv_norm*ex(1:3)
                            cutLineForce(1:3,n)=ex(1:3)
                    END IF
                    ey(1:3) = CrossProduct( ez(1:3), ex(1:3) ) 

                    rotLineForce(1,1:3,n) = ex(1:3)
                    rotLineForce(2,1:3,n) = ey(1:3)
                    rotLineForce(3,1:3,n) = ez(1:3)

                    ! Inverse matrix
                    inv_rotLineForce(1:3,1:3,n) = Transpose( rotLineForce(1:3,1:3,n) )
                    
                    ! Check that rotlineForce is a rotation matrix
                    IF (MatNorm2( MatMul( inv_rotlineForce(:,:,n), rotlineForce(:,:,n) ) - matId ).GT.1.d-6) THEN
                            WRITE(0,'(a,i0,a)') 'Rotation matrix for line force ', n, ' is not a unitary matrix'
                            STOP '< Lecture_input >'
                    END IF
            END DO loop_force1
    ELSE
            nlf=0
    END IF


    !===========================================
    IF ( Count( (/cubic_elasticity, hexagonal_elasticity, &
        anisotropic_elasticity /) ).GE.2 ) THEN
            WRITE(0,'(a)') 'You need to choose between Cubic_Elasticity, &
                &Hexagonal_Elasticity, and Anisotropic_Elasticity'
            STOP '< Lecture_input >'
    END IF
    ! Elastic constants in Voigt notation for cubic symmetry
    IF (cubic_elasticity) THEN
            IF ( (C11.LE.1.d-6).OR.((C11-C12).LE.0.5d-6).OR.(C44.LE.1.d-6) ) THEN
                    WRITE(0,'(a)') 'C11, C11-C12 and C44 have to be strictly &
                            &positive for cubic elasticity'
                    STOP '< Lecture_input >'
            END IF
            CVoigt(1:6,1:6)=0.d0
            CVoigt(1,1)=C11 ; CVoigt(1,2)=C12 ; CVoigt(1,3)=C12
            CVoigt(2,1)=C12 ; CVoigt(2,2)=C11 ; CVoigt(2,3)=C12
            CVoigt(3,1)=C12 ; CVoigt(3,2)=C12 ; CVoigt(3,3)=C11
            CVoigt(4,4)=C44 ; CVoigt(5,5)=C44 ; CVoigt(6,6)=C44
            cubic_elasticity=.FALSE.
            anisotropic_elasticity=.TRUE.
    ELSE IF (hexagonal_elasticity) THEN
            IF ( (C11.LE.1.d-6).OR.(C12.LE.1.d-6).OR.(C13.LE.1.D-6) &
                .OR.(C33.LE.1.D-6).OR.(C44.LE.1.d-6) ) THEN
                    WRITE(0,'(a)') 'C11, C12, C13, C33 and C44 have to be strictly &
                    &positive for hexagonal elasticity'
                    STOP '< Lecture_input >'
            END IF
            CVoigt(1:6,1:6)=0.d0
            CVoigt(1,1)=C11 ; CVoigt(1,2)=C12 ; CVoigt(1,3)=C13
            CVoigt(2,1)=C12 ; CVoigt(2,2)=C11 ; CVoigt(2,3)=C13
            CVoigt(3,1)=C13 ; CVoigt(3,2)=C13 ; CVoigt(3,3)=C33
            CVoigt(4,4)=C44 ; CVoigt(5,5)=C44 ; CVoigt(6,6)=0.5d0*(C11-C12)
            hexagonal_elasticity=.FALSE.
            anisotropic_elasticity=.TRUE.
    END IF
    IF (anisotropic_elasticity) THEN
            IF (MatNorm2(CVoigt).LT.1.d-6) THEN
                    WRITE(0,'(a)') 'Elastic constants CVoigt(1:6,1:6) need to be &
                        &defined for Anisotropic_Elasticity'
                    STOP '< Lecture_input >'
            END IF
            ! Check if elastic constant matrix is symmetric and definite positive 
            CALL Check_CVoigt(CVoigt)
            CALL MatInv(CVoigt, inv_CVoigt)
    END IF

    !===========================================
    IF (strain) THEN
            IF ( (.NOT.anisotropic_elasticity).AND.(MatNorm2(sStrain).GT.1.d-6) ) THEN
                    WRITE(0,'(a)') 'You need to define elastic constants to impose a homogeneous stress'
                    STOP '< Lecture_input >'
            END IF
            sVoigt(:) = tau*ReShape( (/ sStrain(1,1), sStrain(2,2), sStrain(3,3), &
                0.5d0*(sStrain(2,3)+sStrain(3,2)) , &
                0.5d0*(sStrain(1,3)+sStrain(3,1)) , &
                0.5d0*(sStrain(1,2)+sStrain(2,1)) /), (/ 6 /) )
            eVoigt(:) = MatMul( inv_CVoigt, sVoigt )
            eStrain(:,:) = epsi*eStrain(:,:) + ReShape( &
                (/ eVoigt(1), 0.5d0*eVoigt(6), 0.5d0*eVoigt(5), &
                   0.5d0*eVoigt(6), eVoigt(2), 0.5d0*eVoigt(4), &
                   0.5d0*eVoigt(5), 0.5d0*eVoigt(4), eVoigt(3) /), (/3,3/) )
            IF (MatNorm2(eStrain).LE.1.d-6) THEN
                    strain=.false.
            ENDIF
    ELSE
            eStrain(:,:)=0.d0
    ENDIF

    !===========================================
    IF (translate) THEN
            uTranslate(:) = alat*uTranslate(:)
            IF (Sum(uTranslate(:)**2).LE.1.d-12) translate=.FALSE.
    END IF


    !===========================================
    IF (MatNorm2(pImpurity).GT.1.d-6) THEN
            ! A strain associated to a point-defect is defined
            IF (MatNorm2( pImpurity - Transpose(pImpurity) ).GT.1.d-6) THEN
                    WRITE(0,'(a)') 'The elastic dipole of the impurity, pImpurity(:,:), is not symmetric' 
                    STOP '< Lecture_input >'
            END IF
            impurity=.true.
    END IF
    IF (out_Ebinding) THEN
            IF (.NOT.impurity) THEN
                    WRITE(0,'(a)') 'You need to define the elastic dipole of the impurity pImpurity(:,:)'
                    WRITE(0,'(a)') '  to calculate binding energies with out_Ebinding'
            END IF
    END IF


    !===========================================
    IF (slab) THEN
            IF ( (nxSlab.LE.0).OR.(nySlab.LE.0).OR.(nzSlab.LE.0) ) THEN
                    WRITE(0,'(a)') &
                        'Slab: number of slabs in each direction has to be stricly positive'
                    WRITE(0,'(a,i0)') '  nxSlab = ', nxSlab
                    WRITE(0,'(a,i0)') '  nySlab = ', nySlab
                    WRITE(0,'(a,i0)') '  nzSlab = ', nzSlab
                    STOP '< Lecture_input >'
            ELSE IF ( (nxSlab.EQ.1).AND.(nySlab.EQ.1).AND.(nzSlab.EQ.1) ) THEN
                    slab=.false.
            END IF
    END IF



    !===========================================
    IF ( Count( (/outXyz,outOnlyAtoms,outCfg,outGin,outLisa,outSiesta,outLammps,outPoscar/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, OnlyAtoms, Cfg, Gin, Lisa, Siesta, Lammps, and Poscar formats for output structure'
            STOP '< Lecture_input >'
    END IF

END SUBROUTINE Lecture_input
