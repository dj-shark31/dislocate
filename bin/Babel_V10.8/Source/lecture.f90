MODULE babel_lecture

  ! Name of the input file defining the dislocation
  CHARACTER(len=100) :: input_file

  INTEGER, parameter, private :: verbosity_max=2

CONTAINS

  SUBROUTINE Read_Babel(out)

    USE stringModule
    USE babel_data
    USE disloc_elasticity_ani
    USE LineCouple_elasticity_ani
    USE LineForce_elasticity_ani
    USE DDipoleModule
    USE loopModule
    USE fitModule
    USE fitConstraintModule
    USE math
    USE elasticity_ani
    USE rearrange
    USE slabModule
    USE Symmetry3DModule 
    IMPLICIT NONE
    INTEGER, intent(in) :: out

    REAL(kind(0.d0)), dimension(1:6) :: sVoigt
    REAL(kind(0.d0)) :: volume, alat
    INTEGER :: i, n, inp, io

    CHARACTER(len=20) :: keyword

    !===========================================
    ! Initialization
    alat = 1.d0         ! Lattice parameter (A)

    !===========================================
    ! If input is read from "keyboard", write input in a temporary file
    IF (input_file.EQ.'-') THEN
            input_file="input.babel.temp"
            CALL WriteInput(input_file)
    END IF
            
    ! Read input data
    inp=50
    OPEN(file=input_file,unit=inp,action='read',status='old')

    DO
            READ(inp,*,iostat=io) keyword
            IF (io.NE.0) Exit
            BACKSPACE(inp)
            keyword=stringLowerCase(keyword)    ! Put keyword in lower case
            IF (Trim(keyword).EQ."&input") THEN
                    CALL Lecture_input(inp, alat)
            ELSE IF (Trim(keyword).EQ."&dislo") THEN
                    CALL ReadDislo(inp)
            ELSE IF (Trim(keyword).EQ."&ddipole") THEN
                    CALL ReadDDipole(inp)
            ELSE IF (Trim(keyword).EQ."&linecouple") THEN
                    CALL ReadLineCouple(inp)
            ELSE IF (Trim(keyword).EQ."&loops") THEN
                    CALL ReadLoops(inp)
            ELSE IF (Trim(keyword).EQ."&fit") THEN
                    CALL ReadFit(inp)
            ELSE IF (Trim(keyword).EQ."&constraint") THEN
                    CALL ReadConstraint(inp)
            ELSE
                    WRITE(0,'(a)') 'Keyword read in input file: ', Trim(keyword)
                    STOP "Read_Babel"
            END IF
    END DO

    IF (inp.NE.5) CLOSE(inp)
    IF (DEBUG) verbosity=verbosity_debug

    ! Print version information if verbosity is >= 1
    CALL Version(out)

    ! Multiply distances by lattice vectors
    IF (l_dislo) CALL ScaleDislo(alat)
    IF (l_DDipole) CALL ScaleDDipole(alat)
    IF (l_loop) CALL ScaleLoops(alat)
    IF (l_lineCouple) CALL ScaleLineCouple(alat)
    IF (l_constraint) CALL ScaleConstraint(alat)

    ! Make some rearrangements in definitions and calculate auxiliary properties
    ! (rotation matrices)
    IF (l_dislo) CALL RearrangeDislo(out)
    IF (l_DDipole) CALL RearrangeDDipole()
    IF (l_lineCouple) CALL RearrangeLineCouple(out)
    IF (l_dislo.AND.l_lineCouple) CALL Cross_dislo_lineCouple(out)

    !===========================================
    ! Write parameters on output unit
    IF (verbosity.LT.verbosity_max) RETURN
    WRITE(out,*)
    WRITE(out,'(a)') '==========================='
    WRITE(out,*)
    WRITE(out,'(a)') 'Input structure'
    WRITE(out,'(a,i0)') '  number of atoms: ', im
    IF (im.GT.0) WRITE(out,'(a,3g14.6)') '  gravity center located in: ', &
                SUM(xp(1:3,1:im),2)/dble(im)
    IF (.NOT.at_defined) THEN
            WRITE(out,'(a)') '  unit cell vectors not defined' 
    ELSE
            WRITE(out,'(a)') '  unit cell vectors:'
            WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
            WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
            WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
            volume=abs(MatDet(at))
            IF (im.LE.0) THEN
                    WRITE(out,'(a,g20.12)')  '  corresponding volume ', volume
            ELSE
                    WRITE(out,'(2(a,g20.12),a)')  '  corresponding volume ', volume, &
                        '  ( ', volume/dble(im), ' / atom )'
            END IF
    END IF
    WRITE(out,*)

    IF (refXyz.OR.refCfg.OR.refGin.OR.refLisa.OR.refSiesta.OR.refNDM.OR.refLammps.OR.refPoscar) THEN
            WRITE(out,'(a)') 'Reference structure'
            WRITE(out,'(a,i0)') '  number of atoms: ', im0
            IF (im0.GT.0) WRITE(out,'(a,3g14.6)') '  gravity center located in: ', &
                        SUM(xp0(1:3,1:im0),2)/dble(im0)
            IF (.NOT.at_defined) THEN
                    WRITE(out,'(a)') '  unit cell vectors not defined' 
            ELSE
                    WRITE(out,'(a)') '  unit cell vectors:'
                    WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at0(1:3,1)
                    WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at0(1:3,2)
                    WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at0(1:3,3)
                    volume=abs(MatDet(at0))
                    IF (im0.LE.0) THEN
                            WRITE(out,'(a,g20.12)')  '  corresponding volume ', volume
                    ELSE
                            WRITE(out,'(2(a,g20.12),a)')  '  corresponding volume ', volume, &
                                '  ( ', volume/dble(im0), ' / atom )'
                    END IF
            END IF
            WRITE(out,*)
    END IF

    WRITE(out,'(a)') '==========================='
    WRITE(out,*)
    !!$WRITE(out,nml=input)
    !!$WRITE(out,*)
    WRITE(out,'(a)') 'Input parameters'
    WRITE(out,*)
    IF (duplicate) THEN
            WRITE(out,'(a)') ' Duplication of the unit cell'
            WRITE(out,'(2x,i0,a,3g14.6)') lat(1), ' times in direction at(1:3,1): ', &
               at(1:3,1)
            WRITE(out,'(2x,i0,a,3g14.6)') lat(2), ' times in direction at(1:3,2): ', &
               at(1:3,2)
            WRITE(out,'(2x,i0,a,3g14.6)') lat(3), ' times in direction at(1:3,3): ', &
               at(1:3,3)
            WRITE(out,*)
    END IF
    IF (symmetrize) THEN
            WRITE(out,'(a)') ' Symmetrize input structure according to symmetry operations'
            IF (verbosity.GE.verbosity_max+3) CALL PrintSym3Dgroup(symGroup,out)
            WRITE(out,*)
    END IF
    IF (xNoise.GT.0.d0) THEN
            WRITE(out,'(a,g14.6)') ' Noise added to atomic positions: ', xNoise
    END IF
    IF (rotate) THEN
            WRITE(out,'(a)') ' Rotation of the unit cell according to the matrix'
            WRITE(out,'(a,3g14.6,a)') '   | ', rot(1,1:3), ' |'
            WRITE(out,'(a,3g14.6,a)') '   | ', rot(2,1:3), ' |'
            WRITE(out,'(a,3g14.6,a)') '   | ', rot(3,1:3), ' |'
            WRITE(out,*)
    END IF
    IF (translate) THEN
            WRITE(out,'(a)') ' Translation of the unit cell according to the vector'
            WRITE(out,'(a,3g14.6)') '   u(1:3) = ', uTranslate(1:3)
            WRITE(out,*)
    END IF
    IF ( anisotropic_elasticity.OR.isotropic_elasticity ) THEN
            WRITE(out,'(a)') ' Elastic constants in Voigt notation (GPa)'
            CALL Print_CVoigt(CVoigt,out)
            WRITE(out,'(a)') ' Inverse matrix (GPa^-1)'
            CALL Print_CVoigt(inv_CVoigt,out)
            IF ( Isotropic_CVoigt(CVoigt) ) THEN
                    WRITE(out,'(a)') '    elastic constants are isotropic'
                    WRITE(out,'(a,g14.6)') '      lambda = ', CVoigt(1,2)
                    WRITE(out,'(a,g14.6)') '      mu     = ', CVoigt(4,4)
                    WRITE(out,'(a,g14.6)') '      nu     = ', 0.5d0*CVoigt(1,2)/(CVoigt(1,2)+CVoigt(4,4))
            ELSE
                    WRITE(out,'(a)') '    elastic constants are not isotropic'
            END IF
            IF (CVoigt_noise.GT.0.d0) THEN
                    WRITE(out,'(a,g14.6)') '    relative noise to add to elastic&
                        & constants for line defects: ', CVoigt_noise
            END IF
            WRITE(out,*)
    END IF
    IF (l_loop) CALL PrintLoops(out)
    IF (l_dislo) THEN
            WRITE(out,'(a)') ' Dislocations created from anisotropic elasticity theory'
            CALL PrintDislo(out)
    END IF
    IF (l_DDipole) THEN
            WRITE(out,'(a)') ' Dislocation dipoles created from anisotropic elasticity theory'
            CALL PrintDDipole(out)
    END IF
    IF (LineForce) THEN
            WRITE(out,'(a)') ' Line-force created from anisotropic elasticity theory'
            CALL PrintLineForce(out)
    END IF
    IF (l_LineCouple) THEN
            WRITE(out,'(a)') ' Line-force dipole created from anisotropic elasticity theory'
            CALL PrintLineCouple(out)
    END IF
    IF (Strain) THEN
            WRITE(out,'(a)') '  Homogeneous strain applied:'
            DO i=1, 3
                WRITE(out,'(a,3g14.6,a)') '      | ', eStrain(i,1:3), ' |'
            END DO
            WRITE(out,*)
            IF (anisotropic_elasticity) THEN
                    sVoigt = MatMul( CVoigt, (/ eStrain(1,1), eStrain(2,2), estrain(3,3), &
                        eStrain(2,3)+eStrain(3,2), eStrain(1,3)+eStrain(3,1), eStrain(1,2)+eStrain(2,1) /) )
                    WRITE(out,'(a)') '  Corresponding stress (GPa):'
                    WRITE(out,'(a,3g14.6,a)') '      | ', sVoigt(1), sVoigt(6), sVoigt(5), ' |'
                    WRITE(out,'(a,3g14.6,a)') '      | ', sVoigt(6), sVoigt(2), sVoigt(4), ' |'
                    WRITE(out,'(a,3g14.6,a)') '      | ', sVoigt(5), sVoigt(4), sVoigt(3), ' |'
                    WRITE(out,*)
            END IF
    END IF

    IF (xImages) THEN
            WRITE(out,'(a,i0,a)') '  Periodicity in x direction (', &
                nxImages, ' images considered)'
            WRITE(out,'(a,3g14.6)') '    periodicity vector: at(1:3,1) = ', &
                Abs(dble(lat(1)))*at(1:3,1)
    ELSE
            WRITE(out,'(a,i0,a)') '  No periodicity in x direction (', & 
                nxImages, ' image considered)'
    END IF
    IF (yImages) THEN
            WRITE(out,'(a,i0,a)') '  Periodicity in y direction (', &
                nyImages, ' images considered)'
            WRITE(out,'(a,3g14.6)') '    periodicity vector: at(1:3,2) = ', &
                Abs(dble(lat(2)))*at(1:3,2)
    ELSE
            WRITE(out,'(a,i0,a)') '  No periodicity in y direction (', &
                nyImages, ' image considered)'
    END IF
    IF (zImages) THEN
            WRITE(out,'(a,i0,a)') '  Periodicity in z direction (', &
                nzImages, ' images considered)'
            WRITE(out,'(a,3g14.6)') '    periodicity vector: at(1:3,3) = ', &
                Abs(dble(lat(3)))*at(1:3,3)
    ELSE
            WRITE(out,'(a,i0,a)') '  No periodicity in z direction (', &
                nzImages, ' image considered)'
    END IF
    IF (remove_cut) THEN
            WRITE(out,'(a)') '  Atoms necessary to create dislocation or loop cut &
                &will be removed'
    END IF

    IF (slab) THEN
            WRITE(out,'(a)') '  A slab will be extracted from unit cell'
            WRITE(out,'(a,3(i0,a))') '  periodicity vectors will be divided by ' , &
                nxSlab, ', ', nySlab, ', and ', nzSlab, ' respectively'
    END IF
    IF (clipDisplacement) THEN
            WRITE(out,'(a)') '  Apply periodic boundary conditions to &
                &atom displacements'
    ELSE
            WRITE(out,'(a)') '  Periodic boundary conditions to &
                &atom displacements are not applied'
    END IF
    IF (clipAtom) WRITE(out,'(a)') '  Apply periodic boundary conditions to &
                &atom coordinates'
    IF (strain.OR.l_dislo.OR.l_DDipole.OR.l_loop.OR.LineForce.OR.l_LineCouple) THEN
            IF (fixGravity) THEN
                    WRITE(out,'(a)') '  Gravity center is kept fixed when &
                        &creating line-defects and applying homogeneous strain'
            END IF
            IF (initial) THEN
                    WRITE(out,'(a)') '  Keep initial atom coordinates'
            ELSE
                    WRITE(out,'(a)') '  Add displacement field &
                        &to initial atom coordinates'
            END IF
            WRITE(out,*)
    END IF
    IF (max_Euler.EQ.0) THEN
            WRITE(out,'(a)') '  Lagrangian coordinates are used (max_Euler=0)'
    ELSE
            WRITE(out,'(a)') '  Eulerian coordinates are used'
            WRITE(out,'(a,i0,a)')  '    maximal number of iterations: max_Euler=', max_Euler, ')'
            WRITE(out,'(a,g14.6)') '    maximal tolerance for displacement convergency:&
                & delta_Euler = ', delta_Euler
    END IF
    IF (out_Ebinding) THEN
            WRITE(out,'(a)') '  Calculate elastic interaction with impurity'
            WRITE(out,'(a)') '    impurity elastic dipole tensor (eV):'
            DO i=1, 3
                WRITE(out,'(a,3g14.6,a)') '      | ', pImpurity(i,1:3), ' |'
            END DO
            WRITE(out,*)
    END IF
    WRITE(out,'(a)') '==========================='
    WRITE(out,*)

  END SUBROUTINE Read_Babel

END MODULE babel_lecture
