MODULE xyz_module
  ! Structure files in xyz format (compatible with rasmol and glmol)

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter, private :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter, private :: Distance_Zero2=Distance_Zero*Distance_Zero

  PRIVATE :: Find_iTyp
  

CONTAINS

  FUNCTION GetImmXyz(inp) RESULT(imm)
    ! Get maximal number of atoms imm from xyz file connected to unit inp
    
    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    INTEGER :: imm

    imm=0
    READ(inp,*) imm
    REWIND(inp)

  END FUNCTION GetImmXyz

  SUBROUTINE ReadXyz(xp, iTyp, im, at, nTypes, label, inp, nAux, aux)
    ! Read structure in input file connected to inp
    !  Additional columns can be used to read auxiliary properties
    !   nAux: number of auxiliary properties to read
    !   aux(1:nAux,i): corresponding property for atom i

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(out) ::xp
    INTEGER, dimension(:), intent(out) :: iTyp
    INTEGER, intent(out) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: at
    INTEGER, intent(inout) :: nTypes
    CHARACTER(len=5), dimension(:), intent(inout) :: label
    INTEGER, intent(in) :: inp
    INTEGER, intent(in), optional :: nAux
    REAL(kind(0.d0)), dimension(:,:), intent(out), optional :: aux

    INTEGER :: i, io, max_nTypes
    REAL(kind(0.d0)) :: alat
    REAL(kind(0.d0)), dimension(1:3) :: xpTemp
    CHARACTER(len=5) :: symbol
    LOGICAL :: read_aux

    ! Initialization
    max_nTypes=Size(label(:),1)
    at(:,:)=0.d0

    READ(inp,*) im      ! Number of atoms
    READ(inp,*)         ! Title

    ! Maximal number of atoms in simulation box
    IF ( (im.GT.size(xp,2)).OR.(im.GT.size(iTyp,1)) ) THEN
             WRITE(0,'(a,i0)') '  dimension of xp(1:3,:): ', size(xp,2)
             WRITE(0,'(a,i0)') '  dimension of iTyp(:): ', size(iTyp,1)
             WRITE(0,'(a,i0)') '  number of atoms read in XYZ file: ', im
             STOP '< ReadXyz >'
    END IF

    ! Check if auxiliary properties have to be read
    !   and are properly defined
    read_aux=.FALSE.
    IF (Present(nAux)) THEN
            IF (nAux.GT.0) THEN
                    IF (.NOT.Present(aux)) THEN
                            WRITE(0,'(a)') 'Array aux is missing'
                            STOP '< ReadXyz >'
                    END IF
                    IF ( (size(aux,1).LT.nAux).OR.(size(aux,2).LT.im) ) THEN
                            WRITE(0,'(a)') 'Array aux does not have the right size'
                            STOP '< ReadXyz >'
                    END IF
                    read_aux=.TRUE.
            END IF
    END IF

    ! Read atom cartesian coordinates
    IF (read_aux) THEN
            DO i=1, im
               CALL Comment(inp)
               READ(inp,*) symbol, xp(1:3,i), aux(1:nAux,i)
               iTyp(i)=find_iTyp(symbol, nTypes, label, max_nTypes)
            END DO
    ELSE
            DO i=1, im
               CALL Comment(inp)
               READ(inp,*) symbol, xp(1:3,i)
               iTyp(i)=find_iTyp(symbol, nTypes, label, max_nTypes)
            END DO
    END IF

    ! Read periodicity vectors if defined
    CALL Comment(inp)
    io=0
    DO i=1, 3
       CALL Comment(inp)
       READ(inp,*,iostat=io) at(1:3,i)
       IF (io.NE.0) THEN
               at(:,:)=0.d0
               EXIT
       END IF
    END DO
    IF (io.NE.0) THEN
            ! Try to read periodicity vectors in title line (NDM rasmol format)
            REWIND(inp)
            READ(inp,*) ! Number of atoms
            !READ(inp,'(9f12.6)',iostat=io) at(1:3,1), at(1:3,2), at(1:3,3)
            READ(inp,*,iostat=io) at(1:3,1), at(1:3,2), at(1:3,3)
            IF (io.EQ.0) THEN
                    ! There should remain im lines containing atom coordinates
                    DO i=1, im
                       CALL Comment(inp) 
                       READ(inp,*,iostat=io) symbol, xpTemp(1:3)
                       IF (io.NE.0) exit
                    END DO
            END IF
            IF (io.NE.0) at(:,:)=0.d0
            Return
    END IF

    ! Read scaling factor if defined
    CALL Comment(inp)
    READ(inp,*,iostat=io) alat
    IF (io.NE.0) alat=1.D0

    ! Done only in the case where the periodicity vectors are defined
    xp(:,1:im) = alat*xp(:,1:im)
    at(:,:) = alat*at(:,:)

  END SUBROUTINE ReadXyz

  SUBROUTINE WriteXyz(xp, iTyp, im, at, nTypes, label, out_alat, out, &
        mask, nAux_int, aux_int, nAux_real, aux_real, aux_title)
    ! Write structure in output file connected to unit out using Xyz (rasmol) format
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written
    !  Additional columns are used to print auxiliary properties
    !   nAux_int, nAux_real: number of auxiliary properties to print (int and real types)
    !   aux_int(1:nAux_int,i), aux_real(1:nAux_real,i): corresponding property for atom i
    !   aux_title(1:nAux_int): title for auxiliary property with integer type
    !   aux_title(nAux_int+1:nAux_int+nAux_real): title for auxiliary property with real type

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
    INTEGER, intent(in), optional :: nAux_int, nAux_real
    INTEGER, dimension(:,:), optional :: aux_int
    REAL(kind(0.d0)), dimension(:,:), optional :: aux_real
    CHARACTER(len=50), dimension(:), intent(in), optional :: aux_title

    INTEGER :: i, n
    REAL(kind(0.d0)) :: alat
    CHARACTER(len=50) :: out_format
    LOGICAL, dimension(:), allocatable :: local_mask
    LOGICAL :: test_aux_int, test_aux_real, label_defined

    ! Check if atom types have been defined
    label_defined=.FALSE.
    DO n=1, nTypes
       IF (label(n).NE.aChar(n+64)) THEN
               label_defined=.true.
               EXIT
       END IF
    END DO
    IF (.NOT.label_defined) THEN
            WRITE(0,'(a)') 'WARNING < WriteXyz >'
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
    IF ( Present(nAux_int) ) THEN
            IF ( (.NOT.Present(aux_int)) .OR. (.NOT.Present(aux_title)) ) THEN
                    WRITE(0,'(a)') 'All auxiliary properties have to be defined'
                    STOP '< WriteXyz >'
            END IF
            IF (Size(aux_int,1).LT.nAux_int)  THEN
                    WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                    STOP '< WriteXyz >'
            END IF
            test_aux_int=.TRUE.
    ELSE
            test_aux_int=.FALSE.
    END IF
    IF ( Present(nAux_real) ) THEN
            IF ( (.NOT.Present(aux_real)) .OR. (.NOT.Present(aux_title)) ) THEN
                    WRITE(0,'(a)') 'All auxiliary properties have to be defined'
                    STOP '< WriteXyz >'
            END IF
            IF (Size(aux_real,1).LT.nAux_real)  THEN
                    WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                    STOP '< WriteXyz >'
            END IF
            test_aux_real=.TRUE.
    ELSE
            test_aux_real=.FALSE.
    END IF

    ! Parameters used to normalized distances
    IF (out_alat.GT.0.d0) THEN
            alat=out_alat
    ELSE
            alat=1.d0
    END IF


    WRITE(out,'(i0)') Count(local_mask(1:im))        ! Number of atoms

    IF (test_aux_int.AND.test_aux_real) THEN
            WRITE(out_format,'(a,i0,a)') '(a,', nAux_int+nAux_real,'(", ",a))'
            WRITE(out,out_format) 'label, x, y, z', &
                (Trim(AdjustL(aux_title(i))),i=1,nAux_int+nAux_real) ! Title
            WRITE(out_format,'(a,3(i0,a))') &
                '(a,1x,',3,'(g24.16,1x),',nAux_int,'(1x,i0),',nAux_real,'(1x,g24.16))'
            DO i=1, im
               IF (local_mask(i)) WRITE(out,out_format) label(iTyp(i)), xp(1:3,i)/alat, &
                        aux_int(1:nAux_int,i),  aux_real(1:nAux_real,i)
            END DO
    ELSEIF (test_aux_int) THEN
            WRITE(out_format,'(a,i0,a)') '(a,', nAux_int,'(", ",a))'
            WRITE(out,out_format) 'label, x, y, z', &
                (Trim(AdjustL(aux_title(i))),i=1,nAux_int) ! Title
            WRITE(out_format,'(a,2(i0,a))') &
                '(a,1x,',3,'(g24.16,1x),',nAux_int,'(1x,i0))'
            DO i=1, im
               IF (local_mask(i)) WRITE(out,out_format) label(iTyp(i)), xp(1:3,i)/alat, &
                        aux_int(1:nAux_int,i)
            END DO
    ELSEIF (test_aux_real) THEN
            WRITE(out_format,'(a,i0,a)') '(a,', nAux_real,'(", ",a))'
            WRITE(out,out_format) 'label, x, y, z', &
                (Trim(AdjustL(aux_title(i))),i=1,nAux_real) ! Title
            WRITE(out_format,'(a,2(i0,a))') &
                '(a,1x,',3,'(g24.16,1x),',nAux_real,'(1x,g24.16))'
            DO i=1, im
               IF (local_mask(i)) WRITE(out,out_format) label(iTyp(i)), xp(1:3,i)/alat, &
                        aux_real(1:nAux_real,i)
            END DO
    ELSE
            WRITE(out,'(a)') 'label, x, y, z'                ! Title
            WRITE(out_format,'(a,i0,a)') '(a,1x,',3,'(g24.16,1x))'
            DO i=1, im
               IF (local_mask(i)) WRITE(out,out_format) label(iTyp(i)), xp(1:3,i)/alat
            END DO
    END IF

    IF ( SUM( at(1:3,1:3)**2 ).GT.Distance_Zero2 ) THEN 
            WRITE(out,*)
            DO i=1, 3
               WRITE(out,'(3g24.16)') at(1:3,i)/alat ! Basis vector i
            END DO
            WRITE(out,*)
            WRITE(out,'(g24.16)') alat  ! Scaling factor for atom coordinates and basis vectors
    ELSE IF (Abs(alat-1.d0).GT.distance_zero) THEN
            WRITE(out,*)
            DO i=1, 3
               WRITE(out,'(3a)') '? ' ! Basis vector i
            END DO
            WRITE(out,*)
            WRITE(out,'(g24.16)') alat  ! Scaling factor for atom coordinates and basis vectors
    END IF

  END SUBROUTINE WriteXyz

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
      IF (symbol.EQ.label(i)) THEN
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

END MODULE xyz_module
