MODULE drag_lecture

  ! Name of the input file defining the dislocation
  CHARACTER(len=100), save :: input_file

  ! Name of the file where the constraint is stored
  CHARACTER(len=100), save :: constrFile

  ! Reaction coordinates
  REAL(kind(0.d0)), save :: zeta

  ! Number of replica of the simulation box in each direction
  INTEGER, save :: nxSlab, nySlab, nzSlab

  ! Tell if the constraint has to be applied on reduced or cartesian coordinates
  LOGICAL, save :: reduced

  INTEGER, parameter, private :: verbosity_max=2

CONTAINS

  SUBROUTINE Read_Drag(out)

    USE babel_data
    USE structure_module
    USE math
    IMPLICIT NONE
    INTEGER, intent(in) :: out
    REAL(kind(0.d0)) :: alat
    REAL(kind(0.d0)) :: at_norm2, volume0, volume1
    INTEGER :: i

    ! Format of the structure files (input, output, and reference)
    LOGICAL :: inp0Xyz, inp0Cfg, inp0Gin, inp0Lisa, inp0Siesta, inp0NDM, inp0Lammps, inp0Poscar 
    LOGICAL :: inp1Xyz, inp1Cfg, inp1Gin, inp1Lisa, inp1Siesta, inp1NDM, inp1Lammps, inp1Poscar
    ! Input and output structure files
    CHARACTER(len=100) :: inp0File, inp1File

    CHARACTER(len=6) :: fileType


    NAMELIST /input/ &
        inp0File, inp0Xyz, inp0Cfg, inp0Gin, inp0Lisa, inp0Siesta, inp0NDM, inp0Lammps, inp0Poscar, &
        inp1File, inp1Xyz, inp1Cfg, inp1Gin, inp1Lisa, inp1Siesta, inp1NDM, inp1Lammps, inp1Poscar, &
        imm, &
        clipAtom, clipDisplacement, & 
        xNoise, alat, &
        outFile, outXyz, outCfg, outGin, outLisa, outSiesta, outNDM, outLammps, outPoscar, &
        out_alat, &
        constrFile, zeta, reduced, &
        nxSlab, nySlab, nzSlab, &
        nTypes, mass, label, &
        verbosity

    !===========================================
    ! Initialisation
    CALL Random_Seed()

    at_defined=.FALSE.

    fixGravity=.FALSE.
    clipAtom=.FALSE.
    clipDisplacement=.FALSE.
    xNoise=0.d0        ! Noise to add to atomic positions
    alat=1.d0
    at0(:,:)=0.d0; at1(:,:)=0.d0
    inp0Xyz=.FALSE. ; inp0Cfg=.FALSE. ; inp0Gin=.FALSE. ; inp0Lisa=.FALSE. ; inp0Siesta=.FALSE. ; inp0NDM=.FALSE. ; inp0Lammps=.FALSE. ; inp0Poscar=.FALSE.
    inp1Xyz=.FALSE. ; inp1Cfg=.FALSE. ; inp1Gin=.FALSE. ; inp1Lisa=.FALSE. ; inp1Siesta=.FALSE. ; inp1NDM=.FALSE. ; inp1Lammps=.FALSE. ; inp1Poscar=.FALSE.
    outXyz=.FALSE. ; outCfg=.FALSE. ; outGin=.FALSE. ; outLisa=.FALSE. ; outSiesta=.FALSE. ; outNDM=.FALSE. ; outLammps=.FALSE. ; outPoscar=.FALSE.
    inp0File='-' ; inp1File='-' ; outFile='-' 
    constrFile=''
    zeta=0.d0
    imm=0
    verbosity=4
    nxSlab=1
    nySlab=1
    nzSlab=1
    reduced=.FALSE.
    CALL InitStructure()        ! Initialize out_alat, nTypes, mass(:), label(:)

    !===========================================
    ! Read input data
    IF (input_file.EQ.'-') THEN
            READ(5,nml=input)
    ELSE
            OPEN(file=input_file,unit=50,action='read',status='old')
            READ(50,nml=input)
            CLOSE(50)
    END IF

    IF (xNoise.GT.0.d0) xNoise=xNoise*alat

    ! Print version information if verbosity is >= 1
    CALL Version(out)

    IF (program_name.NE."prepareDrag") THEN
            WRITE(0,'(2a)') 'program_name = ', program_name
            STOP '< Read_Drag >'
    END IF

    !===========================================
    IF ( Count( (/inp0Xyz,inp0Cfg,inp0Gin,inp0Lisa,inp0Siesta,inp0NDM,inp0Lammps,inp0Poscar/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin, Lisa, Siesta, NDM, Lammps, and Poscar formats for initial input structure'
            STOP '< Read_Drag >'
    END IF
    IF ( Count( (/inp1Xyz,inp1Cfg,inp1Gin,inp1Lisa,inp1Siesta,inp1NDM,inp1Lammps,inp1Poscar/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin, Lisa, Siesta, NDM, Lammps, and Poscar formats for finaal input structure'
            STOP '< Read_Drag >'
    END IF
    IF ( Count( (/outXyz,outCfg,outGin,outLisa,outSiesta,outNDM,outLammps,outPoscar/) ).GT.1 ) THEN
            WRITE(0,'(a)') 'Choose between Xyz, Cfg, Gin Lisa, Siesta, NDM, Lammps, and Poscar formats for output structure'
            STOP '< Read_Drag >'
    END IF
    !===========================================
    ! Read initial input structure file
    IF (inp0Xyz.OR.inp0Cfg.OR.inp0Gin.OR.inp0Lisa.OR.inp0Siesta.OR.inp0NDM.OR.inp0Lammps.OR.inp0Poscar) THEN
            IF (inp0Xyz)    fileType='xyz'
            IF (inp0Cfg)    fileType='cfg'
            IF (inp0Gin)    fileType='gin'
            IF (inp0Siesta) fileType='siesta'
            IF (inp0Lisa)   fileType='lisa'
            IF (inp0NDM)    fileType='ndm'
            IF (inp0Lammps) fileType='lammps'
            IF (inp0Poscar) fileType='poscar'
            CALL ReadStructure(xp0, iTyp0, im0, imm, at0, inp0File, fileType)
    END IF

    !===========================================
    ! Check if basis vectors have been defined
    IF (SUM( at0(1:3,1:3)**2 ).GT.Distance_Zero2) THEN
            at_defined=.TRUE.
    ELSE
            at_defined=.FALSE.
    END IF
    
    !===========================================
    ! Read final input structure file
    IF (inp1Xyz.OR.inp1Cfg.OR.inp1Gin.OR.inp1Lisa.OR.inp1Siesta.OR.inp1NDM.OR.inp1Lammps.OR.inp1Poscar) THEN
            IF (inp1Xyz)    fileType='xyz'
            IF (inp1Cfg)    fileType='cfg'
            IF (inp1Gin)    fileType='gin'
            IF (inp1Siesta) fileType='siesta'
            IF (inp1Lisa)   fileType='lisa'
            IF (inp1NDM)    fileType='ndm'
            IF (inp1Lammps) fileType='lammps'
            IF (inp1Poscar) fileType='poscar'
            CALL ReadStructure(xp1, iTyp1, im1, imm, at1, inp1File, fileType)
    END IF

    ! Check if basis vectors have been defined
    at_norm2=SUM( at1(1:3,1:3)**2 )
    IF ( (at_norm2.GT.Distance_Zero2).AND. .NOT. at_defined) THEN
            WRITE(0,'(a)') 'Basis vector have been defined for final structure &
                & but not for initial one'
            STOP '< lecture_drag >'
    ELSEIF ( .NOT.(at_norm2.GT.Distance_Zero2).AND.  at_defined) THEN
            WRITE(0,'(a)') 'Basis vector have been defined for initial structure &
                & but not for final one'
            STOP '< lecture_drag >'
    END IF

    ! Check number of atoms
    IF (im0.NE.im1) THEN
            WRITE(0,'(a)') 'Both input structure files do not contain the&
                  & same number of atoms'
            WRITE(0,'(a,i0)') '  im0 = ', im0
            WRITE(0,'(a,i0)') '  im1 = ', im1
            STOP "< Read_Drag >"
    ELSE
            im=im0
    END IF

    ! Check that all atoms are of the same type
    IF (.NOT.Allocated(iTyp0)) THEN
            WRITE(0,'(a)') 'Arry iTyp0 has not been allocated'
            STOP '< Read_Drag >'
    END IF
    IF (.NOT.Allocated(iTyp1)) THEN
            WRITE(0,'(a)') 'Arry iTyp1 has not been allocated'
            STOP '< Read_Drag >'
    END IF
    IF (Allocated(iTyp)) Deallocate(iTyp)
    Allocate(iTyp(1:imm))
    DO i=1, im
       IF (iTyp0(i).nE.iTyp1(i)) THEN
               WRITE(0,'(a,i0,a)') 'Atom ', i, ' is not of the same type in both input files'
               STOP '< Read_Drag >'
       ELSE
               iTyp(1:imm)=iTyp0(1:imm)
       END IF
    END DO
    DEALLOCATE(iTyp0) ; DEALLOCATE(iTyp1)
    
    IF ( (nxSlab.LE.0).OR.(nySlab.LE.0).OR.(nzSlab.LE.0) ) THEN
            WRITE(0,'(a)') &
                'Slab: number of slabs in each direction has to be stricly positive'
            WRITE(0,'(a,i0)') '  nxSlab = ', nxSlab
            WRITE(0,'(a,i0)') '  nySlab = ', nySlab
            WRITE(0,'(a,i0)') '  nzSlab = ', nzSlab
            STOP '< Read_Drag >'
    END IF
    !===========================================
    ! Write parameters on output unit
    IF (verbosity.LT.verbosity_max) RETURN
    WRITE(out,*)
    WRITE(out,'(a)') '==========================='
    WRITE(out,*)
    WRITE(out,'(a)') 'Input corresponding to initial structure'
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
            volume0=abs(MatDet(at0))
            IF (im0.LE.0) THEN
                    WRITE(out,'(a,g20.12)')  '  corresponding volume ', volume0
            ELSE
                    WRITE(out,'(2(a,g20.12),a)')  '  corresponding volume ', volume0, &
                        '  ( ', volume0/dble(im0), ' / atom )'
            END IF
    END IF
    WRITE(out,*)

    WRITE(out,'(a)') 'Input corresponding to final structure'
    WRITE(out,'(a,i0)') '  number of atoms: ', im1
    IF (im1.GT.0) WRITE(out,'(a,3g14.6)') '  gravity center located in: ', &
                SUM(xp1(1:3,1:im1),2)/dble(im1)
    IF (.NOT.at_defined) THEN
            WRITE(out,'(a)') '  unit cell vectors not defined' 
    ELSE
            WRITE(out,'(a)') '  unit cell vectors:'
            WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at1(1:3,1)
            WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at1(1:3,2)
            WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at1(1:3,3)
            volume1=abs(MatDet(at1))
            IF (im1.LE.0) THEN
                    WRITE(out,'(a,g20.12)')  '  corresponding volume ', volume1
            ELSE
                    WRITE(out,'(2(a,g20.12),a)')  '  corresponding volume ', volume1, &
                        '  ( ', volume1/dble(im1), ' / atom )'
            END IF
    END IF
    WRITE(out,*)

    IF (xNoise.GT.0.d0) THEN
            WRITE(out,'(a,g14.6)') ' Noise added to atomic positions: ', xNoise
    END IF

    WRITE(out,'(a,g14.6)') 'Reaction coordinate: zeta = ', zeta

    IF (clipDisplacement) WRITE(out,'(a)') '  Apply periodic boundary conditions to &
                &atom displacements'
    IF (clipAtom) WRITE(out,'(a)') '  Apply periodic boundary conditions to &
                &atom coordinates'
    IF (reduced) THEN
            WRITE(out,'(a)') '  Write constraint for reduced coordinates'
    ELSE
            WRITE(out,'(a)') '  Write constraint for cartesian coordinates'
    END IF

    WRITE(out,'(a)') '==========================='
    WRITE(out,*)

  END SUBROUTINE Read_Drag

END MODULE drag_lecture
