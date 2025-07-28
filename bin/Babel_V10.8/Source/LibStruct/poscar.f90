MODULE poscar_module
  ! Structure files in poscar format (compatible with vasp)

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter, private :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter, private :: Distance_Zero2=Distance_Zero*Distance_Zero

  PRIVATE :: Find_iTyp

  INTERFACE WritePoscar
          MODULE PROCEDURE WritePoscar_sorted
          !MODULE PROCEDURE WritePoscar_unsorted
  END INTERFACE WritePoscar
  

CONTAINS

  FUNCTION GetImmPoscar(inp) RESULT(imm)
    ! Get maximal number of atoms imm from poscar file connected to unit inp
    
    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    INTEGER :: imm

    INTEGER :: i, nAtoms
    CHARACTER(len=256) :: phrase_nAtoms

    imm=0

    READ(inp,*)         ! Title
    READ(inp,*)         ! Scaling constant (lattice parameter)
    DO i=1, 3           ! Periodicity vectors
       READ(inp,*) 
    END DO

    ! Read all atom types and number of atoms for each type
    READ(inp,'(a256)')                  ! Atom types
    READ(inp,'(a256)') phrase_nAtoms    ! Number of atoms
    phrase_nAtoms = AdjustL( phrase_nAtoms )
    imm=0

    ! Loop on atom types
    DO
        ! Exit test
        IF (Trim(phrase_nAtoms).EQ."") EXIT 
       
        ! New atom type
        READ(phrase_nAtoms,*) nAtoms
        IF (nAtoms.LT.1) THEN
                WRITE(0,'(a,i0)') '  number of atoms read in POSCAR file: ', nAtoms
                STOP '< GetImmPoscar >'
        END IF
        imm = imm + nAtoms        ! Last atom index

        ! Next value
        i = Index(phrase_nAtoms," ") + 1
        phrase_nAtoms = AdjustL( phrase_nAtoms(i:) )

    END DO

    REWIND(inp)

  END FUNCTION GetImmPoscar

  SUBROUTINE ReadPoscar(xp, iTyp, im, at, nTypes, label, inp)
    ! Read structure in input file connected to inp

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(out) ::xp
    INTEGER, dimension(:), intent(out) :: iTyp
    INTEGER, intent(out) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: at
    INTEGER, intent(inout) :: nTypes
    CHARACTER(len=5), dimension(:), intent(inout) :: label
    INTEGER, intent(in) :: inp

    INTEGER :: i, max_nTypes, nAtoms
    REAL(kind(0.d0)) :: alat
    REAL(kind(0.d0)), dimension(1:3) :: xc
    CHARACTER(len=5) :: symbol
    CHARACTER(len=256) :: phrase_type, phrase_nAtoms, phrase

    ! Initialization
    max_nTypes=Size(label(:),1)
    alat=1.d0
    at(:,:)=0.d0

    READ(inp,*)         ! Title
    READ(inp,*) alat    ! Scaling constant (lattice parameter)
    DO i=1, 3           ! Periodicity vectors
       READ(inp,*) at(1:3,i)
    END DO
    at(:,:) = alat*at(:,:)

    ! Read all atom types and number of atoms for each type
    READ(inp,'(a256)') phrase_type
    phrase_type = AdjustL( phrase_type )
    READ(inp,'(a256)') phrase_nAtoms
    phrase_nAtoms = AdjustL( phrase_nAtoms )
    im=0

    ! Loop on atom types
    DO
        ! Exit test
        IF ( (Trim(phrase_nAtoms).EQ."").OR.(Trim(phrase_type).EQ."") ) EXIT 
       
        ! New atom type
        READ(phrase_nAtoms,*) nAtoms
        IF (nAtoms.LT.1) THEN
                WRITE(0,'(a,i0)') '  number of atoms read in POSCAR file: ', nAtoms
                STOP '< ReadPoscar >'
        END IF
        i  = im + 1             ! First atom index
        im = im + nAtoms        ! Last atom index
        IF ( (im.GT.size(xp,2)).OR.(im.GT.size(iTyp,1)) ) THEN
                 WRITE(0,'(a,i0)') '  dimension of xp(1:3,:): ', size(xp,2)
                 WRITE(0,'(a,i0)') '  dimension of iTyp(:): ', size(iTyp,1)
                 WRITE(0,'(a,i0)') '  number of atoms read in POSCAR file: ', im
                 STOP '< ReadPosCar >'
        END IF
        symbol = phrase_type(1:Index(phrase_type," ")-1)
        iTyp(i:im)=find_iTyp(symbol, nTypes, label, max_nTypes)

        ! Next value
        i = Index(phrase_type," ") + 1
        phrase_type = AdjustL( phrase_type(i:) )
        i = Index(phrase_nAtoms," ") + 1
        phrase_nAtoms = AdjustL( phrase_nAtoms(i:) )

        !WRITE(6,'(3a,i0,a)') '  type ', Trim(symbol), ' :  ', nAtoms, ' atoms' ! DEBUG

    END DO

    READ(inp,*) phrase
    phrase=AdjustL(phrase)
    IF ( (phrase(1:1).EQ."S").OR.(phrase(1:1).EQ."s") ) THEN
            ! Selective dynamics line
            READ(inp,*) phrase
            phrase=AdjustL(phrase)
    END IF

    ! Read atom cartesian coordinates
    IF ( (phrase(1:1).EQ."C").OR.(phrase(1:1).EQ."c") ) THEN
            ! Carthesian coordinates
            DO i=1, im
               CALL Comment(inp)
               READ(inp,*) xp(1:3,i)
            END DO
            xp(:,1:im) = alat*xp(:,1:im)
    ELSE IF ( (phrase(1:1).EQ."D").OR.(phrase(1:1).EQ."d") ) THEN
            ! Direct coordinates
            DO i=1, im
               CALL Comment(inp)
               READ(inp,*) xc(1:3)
               xp(:,i) = MatMul(at(:,:),xc(:))
            END DO
    ELSE
            WRITE(0,'(a)') 'Error when looking for Carthesian | Direct in POSCAR file'
            WRITE(0,'(2a)') '  read: ', AdjustL(phrase)
            STOP '< ReadPosCar >'
    END IF

    REWIND(inp)

  END SUBROUTINE ReadPoscar

  SUBROUTINE WritePoscar_sorted(xp, iTyp, im, at, nTypes, label, out_alat, out, mask)
    ! Write structure in output file connected to unit out using Poscar (vasp) format
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
    INTEGER, dimension(:), intent(in) :: iTyp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    INTEGER, intent(in) :: nTypes
    CHARACTER(len=5), dimension(:), intent(in) :: label
    REAL(kind(0.d0)), intent(in) :: out_alat
    INTEGER, intent(in) :: out
    LOGICAL, dimension(:), intent(in), optional :: mask

    INTEGER :: i, n, nt, nt_max
    INTEGER, dimension(1:100) :: iTyp0, nTyp0
    REAL(kind(0.d0)) :: alat
    LOGICAL, dimension(:), allocatable :: local_mask
    LOGICAL :: label_defined, atom_type_found
    CHARACTER(len=256) :: phrase_type, phrase_nAtoms

    ! Check if atom types have been defined
    label_defined=.FALSE.
    DO n=1, nTypes
       IF (label(n).NE.aChar(n+64)) THEN
               label_defined=.true.
               EXIT
       END IF
    END DO
    IF (.NOT.label_defined) THEN
            WRITE(0,'(a)') 'WARNING < WritePoscar >'
            WRITE(0,'(a)') 'It seems that atom labels have not been defined'
            DO n=1, Max(nTypes,1)
               WRITE(0,'(a,i0,2a)') '  label for atom type ', n,': ', label(n)
            END DO
            WRITE(0,'(a)') 'You can defined different labels using command label(1)="Fe" in Babel input file'
    END IF

    ! Test input parameters
    IF (Allocated(local_mask)) Deallocate(local_mask)
    Allocate(local_mask(1:im))
    IF (Present(mask)) THEN
            local_mask(1:im)=mask(1:im)
    ELSE
            local_mask = .TRUE.
    END IF

    ! Parameters used to normalized distances
    IF (out_alat.GT.0.d0) THEN
            alat=out_alat
    ELSE
            alat=1.d0
    END IF

    WRITE(out,'(a)') '#######'          ! Title
    WRITE(out,'(g24.16)') alat          ! Scaling factor for atom coordinates and basis vectors
    DO i=1, 3
       WRITE(out,'(3g24.16)') at(1:3,i)/alat ! Basis vector i
    END DO

    ! Count number of atoms for each type
    iTyp0(:)=0
    nTyp0(:)=0
    nt_max=0
    DO i=1, im
       atom_type_found=.false.
       IF (local_mask(i)) THEN
               DO nt=1, nt_max
                  IF (iTyp(i).EQ.iTyp0(nt)) THEN
                          ! One more atom for this type
                          nTyp0(nt)=nTyp0(nt)+1
                          atom_type_found=.true.
                  END IF
               END DO
               IF (.NOT.atom_type_found) THEN
                       ! New type of atom
                       nt_max = nt_max + 1
                       IF (nt_max.GT.100) STOP '< WritePoscar >: increase nt_max'
                       nTyp0(nt_max) = 1
                       iTyp0(nt_max) = iTyp(i)
               END IF
       END IF
    END DO

    ! Print number of atoms for each type
    phrase_type=''
    phrase_nAtoms=''
    DO nt=1, nt_max
       WRITE(phrase_type,'(3a)') Trim(phrase_type), " ", Trim(label(iTyp0(nt)))
       WRITE(phrase_nAtoms,'(2a,i0)') Trim(phrase_nAtoms), " ", nTyp0(nt)
    END DO
    ! Write atom types and corresponding number of atoms
    WRITE(out,'(a)') Trim(phrase_type)
    WRITE(out,'(a)') Trim(phrase_nAtoms)

    ! Write atom coordinates
    WRITE(out,'(a)') 'Cartesian'
    DO nt=1, nt_max
       DO i=1, im
          IF ( (local_mask(i)).AND.(iTyp(i).EQ.iTyp0(nt)) ) &
                  WRITE(out,'(3(g24.14,1x))') xp(1:3,i)/alat
       END DO
    END DO

  END SUBROUTINE WritePoscar_sorted

  SUBROUTINE WritePoscar_unsorted(xp, iTyp, im, at, nTypes, label, out_alat, out, mask)
    ! Write structure in output file connected to unit out using Poscar (vasp) format
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written
    !  atoms are not gathered by atom types, keeping initial order

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
    INTEGER, dimension(:), intent(in) :: iTyp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    INTEGER, intent(in) :: nTypes
    CHARACTER(len=5), dimension(:), intent(in) :: label
    REAL(kind(0.d0)), intent(in) :: out_alat
    INTEGER, intent(in) :: out
    LOGICAL, dimension(:), intent(in), optional :: mask

    INTEGER :: i, n, iTyp0, nTyp0
    REAL(kind(0.d0)) :: alat
    LOGICAL, dimension(:), allocatable :: local_mask
    LOGICAL :: label_defined
    CHARACTER(len=256) :: phrase_type, phrase_nAtoms

    ! Check if atom types have been defined
    label_defined=.FALSE.
    DO n=1, nTypes
       IF (label(n).NE.aChar(n+64)) THEN
               label_defined=.true.
               EXIT
       END IF
    END DO
    IF (.NOT.label_defined) THEN
            WRITE(0,'(a)') 'WARNING < WritePoscar >'
            WRITE(0,'(a)') 'It seems that atom labels have not been defined'
            DO n=1, Max(nTypes,1)
               WRITE(0,'(a,i0,2a)') '  label for atom type ', n,': ', label(n)
            END DO
            WRITE(0,'(a)') 'You can defined different labels using command label(1)="Fe" in Babel input file'
    END IF

    ! Test input parameters
    IF (Allocated(local_mask)) Deallocate(local_mask)
    Allocate(local_mask(1:im))
    IF (Present(mask)) THEN
            local_mask(1:im)=mask(1:im)
    ELSE
            local_mask = .TRUE.
    END IF

    ! Parameters used to normalized distances
    IF (out_alat.GT.0.d0) THEN
            alat=out_alat
    ELSE
            alat=1.d0
    END IF

    WRITE(out,'(a)') '#######'          ! Title
    WRITE(out,'(g24.16)') alat          ! Scaling factor for atom coordinates and basis vectors
    DO i=1, 3
       WRITE(out,'(3g24.16)') at(1:3,i)/alat ! Basis vector i
    END DO

    ! Count number of atoms for each type
    phrase_type=''
    phrase_nAtoms=''
    iTyp0=0
    nTyp0=0
    DO i=1, im
       IF (local_mask(i)) THEN
               IF (iTyp(i).EQ.iTyp0) THEN
                       ! One more atom for this type
                       nTyp0=nTyp0+1
               ELSE
                       ! New type of atom
                       IF (nTyp0.NE.0) THEN
                               ! Record previous atom types and corresponding number of atoms
                               WRITE(phrase_type,'(3a)') Trim(phrase_type), " ", Trim(label(iTyp0))
                               WRITE(phrase_nAtoms,'(2a,i0)') Trim(phrase_nAtoms), " ", nTyp0
                               !WRITE(6,'(3a,i0)') 'atom type: ', Trim(label(iTyp0)), &  ! DEBUG
                                       !'  - number of atoms: ', nTyp0                   ! DEBUG
                       END IF
                       iTyp0=iTyp(i)
                       nTyp0=1
               END IF
       END IF
    END DO
    WRITE(phrase_type,'(3a)') Trim(phrase_type), " ", Trim(label(iTyp0))
    WRITE(phrase_nAtoms,'(2a,i0)') Trim(phrase_nAtoms), " ", nTyp0
    !WRITE(6,'(3a,i0)') 'atom type: ', Trim(label(iTyp0)), '  - number of atoms: ', nTyp0 ! DEBUG
    ! Write atom types and corresponding number of atoms
    WRITE(out,'(a)') Trim(phrase_type)
    WRITE(out,'(a)') Trim(phrase_nAtoms)

    ! Write atom coordinates
    WRITE(out,'(a)') 'Cartesian'
    DO i=1, im
       IF (local_mask(i)) WRITE(out,'(3(g24.14,1x))') xp(1:3,i)/alat
    END DO

  END SUBROUTINE WritePoscar_unsorted

  FUNCTION Find_iTyp(symbol, nTypes, label, max_nTypes) RESULT(iAtom)

    IMPLICIT NONE
    CHARACTER(len=5), intent(in) :: symbol
    INTEGER, intent(inout) :: nTypes
    CHARACTER(len=5), dimension(:), intent(inout) :: label
    INTEGER, intent(in) :: max_nTypes
    INTEGER :: iAtom

    INTEGER :: i

    ! Look for symbol in table label
    DO i=1, nTypes
      IF (Trim(symbol).EQ.Trim(label(i))) THEN
              iAtom=i
              RETURN
      END IF
    END DO
    
    ! Program does not find symbol in table label
    !  => we need to define a new type
    nTypes=nTypes+1
    IF (nTypes.GT.max_nTypes) THEN
            WRITE(0,'(a)') 'Too many atom types found'
            WRITE(0,'(a)') 'You need to increase parameter max_nTypes in structure.f90 and to recompile'
            WRITE(0,'(a,i0)') '  current value: ', max_nTypes
            STOP '< Find_iTyp >'
    END IF
    label(nTypes)=symbol
    iAtom=nTypes

  END FUNCTION Find_iTyp

END MODULE poscar_module
