MODULE Structure_Module

  USE MathStruct
  USE Transformation
  IMPLICIT NONE

  SAVE

  ! Parameter used to normalized distance with xyz or cfg output structure file 
  REAL(kind(0.d0)) :: out_alat

  ! Number of different types and corresponding labels and masses
  INTEGER, parameter :: max_nTypes=360
  INTEGER :: nTypes
  CHARACTER(len=5), dimension(1:max_nTypes) :: label
  REAL(kind(0.d0)), dimension(1:max_nTypes) :: mass

CONTAINS

  SUBROUTINE InitStructure()

     IMPLICIT NONE
     INTEGER :: i

     out_alat=1.d0
     nTypes=0
     mass(:)=0.d0
     DO i=1, max_nTypes
       label(i)=aChar(Modulo(i-1,26)+65)
     END DO


  END SUBROUTINE InitStructure

  SUBROUTINE ReadStructure(xp, iTyp, im, imm, at, inpFile, fileType)

     USE cfg_module
     USE xyz_module
     USE gin_module
     USE siesta_module
     USE lisa_module
     USE ndm_module
     USE lammps_module
     USE poscar_module

     IMPLICIT NONE
     ! Atom real coordinates: xp(1:3,:)
     REAL(kind(0.d0)), dimension(:,:), allocatable, intent(out) :: xp
     ! Atom type (input and reference structures)
     INTEGER, dimension(:), allocatable, intent(out) :: iTyp
     ! Number of atoms
     INTEGER, intent(out) :: im
     ! Maximal number of atoms (used to allocate tables)
     INTEGER, intent(inout) :: imm
     ! Lattice vector coordinates (A): at(1:3,1), at(2:3,2), ...
     REAL(kind(0.d0)), dimension(1:3, 1:3), intent(out) :: at
     ! File where to read structure
     CHARACTER(len=*), intent(in) :: inpFile
     ! fileType="cfg", "xyz", "siesta", "gin", "lisa"
     CHARACTER(len=*), intent(in) :: fileType

     INTEGER :: inpUnit


     ! Determine where to read input file
     IF (inpFile.EQ.'-') THEN    
             inpUnit=5           !  from keyboard
     ELSE
             inpUnit=51          !  from file
             IF (Trim(fileType).EQ."ndm") THEN
                     OPEN(file=inpFile,unit=inpUnit,form='unformatted',action='read',status='old')
             ELSE
                     OPEN(file=inpFile,unit=inpUnit,action='read',status='old')
             END IF
     END IF

     ! Get maximal number of atoms if not already defined
     IF (imm.EQ.0) THEN
             IF (inpUnit.EQ.5) THEN
                     WRITE(0,'(a)') &
                         'You need to define the maximal number of atoms, imm, in the input file'
                     WRITE(0,'(a)') &
                         'if you want to read structure file from the keyboard, ie with inpFile="-"'
                     STOP '< ReadStructure >'
             END IF
             SELECT CASE(Trim(fileType))
             CASE("xyz")
                     imm=GetImmXyz(inpUnit)
             CASE("cfg")
                     imm=GetImmCfg(inpUnit)
             CASE("gin")
                     imm=GetImmGin(inpUnit)
             CASE("lisa")
                     imm=GetImmLisa(inpUnit)
             CASE("siesta")
                     imm=GetImmSiesta(inpUnit)
             CASE("ndm")
                     imm=GetImmNDM(inpUnit)
             CASE("lammps")
                     imm=GetImmLammps(inpUnit)
             CASE("poscar")
                     imm=GetImmPoscar(inpUnit)
             CASE DEFAULT
                     WRITE(0,'(2a)') "Unknown format for structure file: ", fileType
                     STOP '< ReadStructure >'
             END SELECT
     END IF

     ! Table allocation
     IF (Allocated(xp)) Deallocate(xp)
     ALLOCATE(xp(1:3,1:imm))
     IF (Allocated(ityp)) Deallocate(ityp)
     ALLOCATE(iTyp(1:imm))

    ! Read input structure file 
     SELECT CASE(Trim(fileType))
     CASE("xyz")
             CALL ReadXyz(xp, iTyp, im, at, nTypes, label, inpUnit)
     CASE("cfg")
             CALL ReadCfg(xp, iTyp, im, at, nTypes, mass, label, inpUnit)
     CASE("gin")
             CALL ReadGin(xp, iTyp, im, at, nTypes, inpUnit)
     CASE("lisa")
             CALL ReadLisa(xp, iTyp, im, at, nTypes, label, inpUnit)
     CASE("siesta")
             CALL ReadSiesta(xp, iTyp, im, at, nTypes, inpUnit)
     CASE("ndm")
             CALL ReadNDM(xp, iTyp, im, at, nTypes, inpUnit)
     CASE("lammps")
             CALL ReadLammps(xp, iTyp, im, at, nTypes, inpUnit)
     CASE("poscar")
             CALL ReadPoscar(xp, iTyp, im, at, nTypes, label, inpUnit)
     CASE DEFAULT
             WRITE(0,'(2a)') "Unknown format for structure file: ", fileType
             STOP '< ReadStructure >'
     END SELECT
            
     ! Close input unit
     IF (inpUnit.NE.5) Close(inpUnit)

  END SUBROUTINE ReadStructure

  SUBROUTINE WriteStructure(xp, iTyp, im, at, outFile, fileType, &
        mask, nAux_int, aux_int, nAux_real, aux_real, aux_title)

     ! Write structure in output file outFile format given in fileType
     !  If mask is given, only atoms for which mask(i)=.TRUE. are written
     !  Additional columns are used to print auxiliary properties
     !   nAux_int, nAux_real: number of auxiliary properties to print (int and real types)
     !   aux_int(1:nAux_int), aux_real(1:nAux_real): corresponding property for atom i
     !   aux_title(1:nAux_int): title for auxiliary property with integer type
     !   aux_title(nAux_int+1:nAux_int+nAux_real): title for auxiliary property with real type

     USE cfg_module
     USE xyz_module
     USE gin_module
     USE siesta_module
     USE lisa_module
     USE ndm_module
     USE onlyAtoms_module
     USE lammps_module
     USE poscar_module

     IMPLICIT NONE
     REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
     INTEGER, dimension(:), intent(in) :: iTyp
     INTEGER, intent(in) :: im
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
     ! File where to write structure
     CHARACTER(len=*), intent(in) :: outFile
     ! fileType="cfg", "xyz", "siesta", "gin", "lisa", "onlyAtoms"
     CHARACTER(len=*), intent(in) :: fileType
     LOGICAL, dimension(:), intent(in), optional :: mask
     INTEGER, intent(in), optional :: nAux_int, nAux_real
     INTEGER, dimension(:,:), intent(in), optional :: aux_int
     REAL(kind(0.d0)), dimension(:,:), intent(in), optional :: aux_real
     CHARACTER(len=50), dimension(:), intent(in), optional :: aux_title

     INTEGER :: outUnit
     LOGICAL :: test_aux_real, test_aux_int

     ! Determine where to write output structure (file or screen)
     IF (outFile.EQ.'-') THEN
             outUnit=6
     ELSE
             IF (Trim(fileType).EQ.'ndm') THEN
                     OPEN(file=outFile,unit=61,action='write',status='unknown',form='unformatted')
             ELSE
                     OPEN(file=outFile,unit=61,action='write',status='unknown')
             END IF
             outUnit=61
     END IF

     ! Test input parameters
     IF ( Present(nAux_int) ) THEN
             IF ( (.NOT.Present(aux_int)) .OR. (.NOT.Present(aux_title)) ) THEN
                     WRITE(0,'(a)') 'All auxiliary properties have to be defined'
                     STOP '< WriteStructure > 1' 
             END IF
             IF (Size(aux_int,1).LT.nAux_int)  THEN
                     WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                     STOP '< WriteStructure > 2'
             END IF
             IF (nAux_int.GT.0) THEN
                     test_aux_int=.TRUE.
             ELSE
                     test_aux_int=.FALSE.
             END IF
     ELSE
             test_aux_int=.FALSE.
     END IF
     IF ( Present(nAux_real) ) THEN
             IF ( (.NOT.Present(aux_real)) .OR. (.NOT.Present(aux_title)) ) THEN
                     WRITE(0,'(a)') 'All auxiliary properties have to be defined'
                     STOP '< WriteStructure > 3'
             END IF
             IF (Size(aux_real,1).LT.nAux_real)  THEN
                     WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                     STOP '< WriteStructure > 4'
             END IF
             IF (nAux_real.GT.0) THEN
                     test_aux_real=.TRUE.
             ELSE
                     test_aux_real=.FALSE.
             END IF
     ELSE
             test_aux_real=.FALSE.
     END IF
     IF (test_aux_int.AND.test_aux_real) THEN
             IF (Size(aux_title,1).LT.nAux_real+nAux_int)  THEN
                     WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                     STOP '< WriteStructure > 5'
             END IF
     END IF

    ! Write output structure file 
     SELECT CASE(Trim(fileType))
     CASE("xyz")
             IF (Present(mask)) THEN
                     IF ( (test_aux_int).AND.(test_aux_real) ) THEN
                             CALL WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, outUnit, mask=mask(:), &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_int) THEN
                             CALL WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, outUnit, mask=mask(:), &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_real) THEN
                             CALL WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, outUnit, mask=mask(:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE
                             CALL WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, outUnit, mask=mask(:))
                     END IF
             ELSE
                     IF ( (test_aux_int).AND.(test_aux_real) ) THEN
                             CALL WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, outUnit, &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_int) THEN
                             CALL WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, outUnit, &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_real) THEN
                             CALL WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, outUnit, &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE
                             CALL WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, outUnit)
                     END IF
             END IF

     CASE("onlyAtoms")
             IF (Present(mask)) THEN
                     IF ( (test_aux_int).AND.(test_aux_real) ) THEN
                             CALL WriteOnlyAtoms(xp, im, out_alat, outUnit, mask=mask(:), &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_int) THEN
                             CALL WriteOnlyAtoms(xp, im, out_alat, outUnit, mask=mask(:), &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_real) THEN
                             CALL WriteOnlyAtoms(xp, im, out_alat, outUnit, mask=mask(:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE
                             CALL WriteOnlyAtoms(xp, im, out_alat, outUnit, mask=mask(:))
                     END IF
             ELSE
                     IF ( (test_aux_int).AND.(test_aux_real) ) THEN
                             CALL WriteOnlyAtoms(xp, im, out_alat, outUnit, &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_int) THEN
                             CALL WriteOnlyAtoms(xp, im, out_alat, outUnit, &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_real) THEN
                             CALL WriteOnlyAtoms(xp, im, out_alat, outUnit, &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE
                             CALL WriteOnlyAtoms(xp, im, out_alat, outUnit)
                     END IF
             END IF

     CASE("cfg")
             IF (Present(mask)) THEN
                     IF ( (test_aux_int).AND.(test_aux_real) ) THEN
                             CALL WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, outUnit, mask=mask(:), &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_int) THEN
                             CALL WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, outUnit, mask=mask(:), &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_real) THEN
                             CALL WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, outUnit, mask=mask(:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE
                             CALL WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, outUnit, mask=mask(:))
                     END IF
             ELSE
                     IF ( (test_aux_int).AND.(test_aux_real) ) THEN
                             CALL WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, outUnit, &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_int) THEN
                             CALL WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, outUnit, &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_real) THEN
                             CALL WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, outUnit, &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE
                             CALL WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, outUnit)
                     END IF
             END IF
             
     CASE("gin")
             IF (Present(mask)) THEN
                     CALL WriteGin(xp, iTyp, im, at, outUnit, mask)
             ELSE
                     CALL WriteGin(xp, iTyp, im, at, outUnit)
             END IF
             
     CASE("lisa")
             IF (Present(mask)) THEN
                     CALL WriteLisa(xp, iTyp, im, at, label, outUnit, mask)
             ELSE
                     CALL WriteLisa(xp, iTyp, im, at, label, outUnit)
             END IF
             
     CASE("siesta")
             IF (Present(mask)) THEN
                     CALL WriteSiesta(xp, iTyp, im, at, outUnit, mask)
             ELSE
                     CALL WriteSiesta(xp, iTyp, im, at, outUnit)
             END IF
             
     CASE("ndm")
             IF (Present(mask)) THEN
                     CALL WriteNDM(xp, iTyp, im, at, outUnit, mask)
             ELSE
                     CALL WriteNDM(xp, iTyp, im, at, outUnit)
             END IF

     CASE("lammps")
             IF (Present(mask)) THEN
                     IF ( (test_aux_int).AND.(test_aux_real) ) THEN
                             CALL WriteLammps(xp, iTyp, im, at, outUnit, mask=mask(:), &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_int) THEN
                             CALL WriteLammps(xp, iTyp, im, at, outUnit, mask=mask(:), &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_real) THEN
                             CALL WriteLammps(xp, iTyp, im, at, outUnit, mask=mask(:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE
                             CALL WriteLammps(xp, iTyp, im, at, outUnit, mask=mask(:))
                     END IF
             ELSE
                     IF ( (test_aux_int).AND.(test_aux_real) ) THEN
                             CALL WriteLammps(xp, iTyp, im, at, outUnit, &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_int) THEN
                             CALL WriteLammps(xp, iTyp, im, at, outUnit, &
                                   nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                   aux_title=aux_title(:))
                     ELSE IF (test_aux_real) THEN
                             CALL WriteLammps(xp, iTyp, im, at, outUnit, &
                                   nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                   aux_title=aux_title(:))
                     ELSE
                             CALL WriteLammps(xp, iTyp, im, at, outUnit)
                     END IF
             END IF
             
     CASE("poscar")
             IF (Present(mask)) THEN
                                     
                     CALL WritePoscar(xp, iTyp, im, at, nTypes, label, out_alat, outUnit, mask)
             ELSE
                     CALL WritePoscar(xp, iTyp, im, at, nTypes, label, out_alat, outUnit)
             END IF

     CASE DEFAULT
             WRITE(0,'(2a)') "Unknown format for structure file: ", fileType
             STOP '< WriteStructure >'
     END SELECT

     ! Close output unit
     IF (outUnit.NE.6) Close(outUnit)


  END SUBROUTINE WriteStructure

END MODULE Structure_Module
