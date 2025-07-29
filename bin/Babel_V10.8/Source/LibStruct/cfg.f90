MODULE cfg_module

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter, private :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter, private :: Distance_Zero2=Distance_Zero*Distance_Zero
  
CONTAINS
  
  FUNCTION GetImmCfg(inp) RESULT(imm)
    ! Get maximal number of atoms imm from cfg file connected to unit inp

    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    INTEGER :: imm

    CHARACTER(len=50) :: phrase
    INTEGER :: pos, io

    imm=0
     ! ===== Preamble =========================
     DO 
        CALL comment(inp)
        READ(inp,'(a50)') phrase

        ! Number of atoms in unit cell
        pos=Index(phrase,'Number of particles')
        IF (pos.NE.0) THEN
                pos=Index(phrase,'=')
                phrase=phrase(pos+1:Len_Trim(phrase))
                READ(phrase,*,iostat=io) imm
                IF (io.NE.0) THEN
                        WRITE(0,'(a)') 'Error when reading in CFG number of particles'
                        WRITE(0,'(2a)') 'Line read in cfg file: ', phrase
                        STOP '< GetImmCfg >'
                END IF
                cycle
        END IF
         
        ! Scaling factor (A)
        pos=Index(phrase,'A =')
        IF (pos.NE.0) cycle

        ! Lattice vector coordinates
        pos=Index(phrase,'H0(')
        IF (pos.NE.0) cycle

        pos=Index(phrase,'.NO_VELOCITY.')
        IF (pos.NE.0) cycle

        ! Number of auxiliary properties
        pos=Index(phrase,'entry_count')
        IF (pos.NE.0) Cycle

        pos=Index(phrase,'auxiliary[')
        IF (pos.NE.0) cycle
        
        BACKSPACE(inp)
        exit
     END DO

     ! Check number of atoms
     IF (imm.EQ.0) THEN
             WRITE(0,*) 'Number of atoms in unit cell is 0'
             WRITE(0,'(a)') 'Line "Number of particles =" should be missing in CFG file'
             STOP '< GetImmCfg >'
     END IF

     REWIND(inp)

  END FUNCTION GetImmCfg

  SUBROUTINE ReadCfg(xp, iTyp, im, at, nTypes, mass, label, inp, nAux, aux)
     ! Read configuration in file *.cfg connected to unit inp

     !USE Math
     IMPLICIT NONE
     REAL(kind(0.d0)), dimension(:,:), intent(out) ::xp
     INTEGER, dimension(:), intent(out) :: iTyp
     INTEGER, intent(out) :: im
     REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: at
     INTEGER, intent(in) :: inp
     INTEGER, intent(out) :: nTypes
     REAL(kind(0.d0)), dimension(:), intent(out) :: mass
     CHARACTER(len=5), dimension(:), intent(out) :: label
     INTEGER, intent(in), optional :: nAux
     REAL(kind(0.d0)), dimension(:,:), intent(out), optional :: aux

     CHARACTER(len=500) :: phrase
     CHARACTER(len=1) :: indice
     INTEGER :: pos, i, j, io, nEntry_Count
     REAL(kind(0.d0)) :: alat
     REAL(kind(0.d0)), dimension(:,:), allocatable :: xc
     

     ! Default values
     im=0           ! Number of atoms in simulation box
     alat=1.d0          ! Scaling factor (A)
     at(1:3,1:3)=0.d0   ! Lattice vector coordinates

     ! Read global parameters
     ! ===== Preamble =========================
     DO 
        CALL comment(inp)
        READ(inp,'(a50)') phrase

        ! Number of atoms in unit cell
        pos=Index(phrase,'Number of particles')
        IF (pos.NE.0) THEN
                pos=Index(phrase,'=')
                phrase=phrase(pos+1:Len_Trim(phrase))
                READ(phrase,*,iostat=io) im
                IF (io.NE.0) THEN
                        WRITE(0,'(a)') 'Error when reading in CFG number of particles'
                        WRITE(0,'(2a)') 'Line read in cfg file: ', phrase
                        STOP '< ReadCfg >'
                END IF
                cycle
        END IF
         
        ! Scaling factor (A)
        pos=Index(phrase,'A =')
        IF (pos.NE.0) THEN
                pos=Index(phrase,'=')
                phrase=phrase(pos+1:Len_Trim(phrase))
                READ(phrase,*,iostat=io) alat
                IF (io.NE.0) THEN
                        WRITE(0,'(a)') 'Error when reading in CFG file scaling factor A'
                        WRITE(0,'(2a)') 'Line read in cfg file: ', phrase
                        STOP '< ReadCfg >'
                END IF
                cycle
        END IF

        ! Lattice vector coordinates
        pos=Index(phrase,'H0(')
        IF (pos.NE.0) THEN
                indice=phrase(pos+3:pos+3)
                READ(indice,*) i
                indice=phrase(pos+5:pos+5)
                READ(indice,*) j
                pos=Index(phrase,'=')
                phrase=phrase(pos+1:Len_Trim(phrase))
                READ(phrase,*,iostat=io) at(j,i)
                IF (io.NE.0) THEN
                        WRITE(0,'(a,2(i0,a))') 'Error when reading in CFG file Lattice vector coordinate H0(', i,',',j,')'
                        WRITE(0,'(2a)') 'Line read in cfg file: ', phrase
                        STOP '< ReadCfg >'
                END IF
                cycle
        END IF

        pos=Index(phrase,'.NO_VELOCITY.')
        IF (pos.NE.0) cycle

        ! Number of auxiliary properties
        pos=Index(phrase,'entry_count')
        IF (pos.NE.0) THEN
                pos=Index(phrase,'=')
                phrase=phrase(pos+1:Len_Trim(phrase))
                READ(phrase,*,iostat=io) nEntry_count
                Cycle
        END IF

        pos=Index(phrase,'auxiliary[')
        IF (pos.NE.0) cycle
        
        BACKSPACE(inp)
        exit
     END DO
     ! ===== Preamble =========================
     
     ! Check lattice parameter
     IF (SUM( at(1:3,1:3)**2 ).LE.Distance_Zero2) THEN
             WRITE(0,'(a)') 'Program does not manage to read lattice parameters in input cfg file'
             STOP '< ReadCfg >'
     END IF

     ! Check number of atoms
     IF (im.EQ.0) THEN
             WRITE(0,*) 'Number of atoms in unit cell is 0'
             WRITE(0,'(a)') 'Line "Number of particles =" should be missing in CFG file'
             STOP '< ReadCfg >'
     END IF

     ! Maximal number of atoms in simulation box (used to allocate tables)
    IF ( (im.GT.size(xp,2)).OR.(im.GT.size(iTyp,1)) ) THEN
             WRITE(0,'(a,i0)') '  dimension of xp(1:3,:): ', size(xp,2)
             WRITE(0,'(a,i0)') '  dimension of iTyp(:): ', size(iTyp,1)
             WRITE(0,'(a,i0)') '  number of atoms read in CFG file: ', im
             STOP '< ReadCfg >'
     END IF

     ! Check number of auxiliary properties
     IF (Present(nAux)) THEN
             IF (nAux.GT.nEntry_Count-3) THEN
                     WRITE(0,'(a)') 'Read in cfg file:'
                     WRITE(0,'(a,i0)')'  entry_count = ', nEntry_Count
                     WRITE(0,'(a,i0,a)') 'Try to read ', nAux, ' auxiliary properties'
                     STOP '< ReadCfg >'
             END IF
     END IF

     ! Table allocation
     IF (Allocated(xc)) Deallocate(xc)
     ALLOCATE(xc(1:3,1:im))

     ! Lattice vectors
     at(1:3,1:3) = alat*at(1:3,1:3)

     ! Initialization: only 1 atom type
     nTypes=1
     ityp(1:im)=1  
     READ(inp,'(a500)',iostat=io) phrase
     IF (io.EQ.0) READ(phrase,*,iostat=io) mass(nTypes)
     IF (io.NE.0) THEN
             WRITE(0,'(a)') 'Error when reading mass in CFG file'
             WRITE(0,'(2a)') 'Line read in cfg file: ', Trim(phrase)
             STOP '< ReadCfg >'
     END IF
     READ(inp,'(a500)',iostat=io) phrase
     IF (io.EQ.0) READ(phrase,*,iostat=io) label(nTypes)
     IF (io.NE.0) THEN
             WRITE(0,'(a)') 'Error when reading atom type in CFG file'
             WRITE(0,'(2a)') 'Line read in cfg file: ', Trim(phrase)
             STOP '< ReadCfg >'
     END IF

     i=0
     DO 
        i=i+1   ! Increment number of atoms
        IF (i.GT.im) EXIT       ! Maximal number of atoms has been reached
        READ(inp,'(a500)',iostat=io) phrase
        IF (io.EQ.0) THEN
                IF (Present(nAux)) THEN
                        READ(phrase,*,iostat=io) xc(1:3,i), aux(1:nAux,i)
                ELSE
                        READ(phrase,*,iostat=io) xc(1:3,i)
                END IF
                IF (io.NE.0) THEN
                        ! Try to see if we do not have a new atom type here
                        nTypes = nTypes + 1
                        !   - read mass of new type
                        READ(phrase,*,iostat=io) mass(nTypes)
                        IF (io.NE.0) THEN
                                WRITE(0,'(a)') 'Error when reading mass for new type in CFG file'
                                WRITE(0,'(a,i0)') '  actual type: ', nTypes
                                WRITE(0,'(2a)') 'Line read in cfg file: ', Trim(phrase)
                                STOP '< ReadCfg >'
                        END IF
                        !   - read label of new type
                        READ(inp,'(a500)',iostat=io) phrase
                        IF (io.EQ.0) READ(phrase,*,iostat=io) label(nTypes)
                        IF (io.NE.0) THEN
                                WRITE(0,'(a)') 'Error when reading atom label for new type in CFG file'
                                WRITE(0,'(a,i0)') '  actual type: ', nTypes
                                WRITE(0,'(a,f0.3)') '  mass read just before : ', mass(nTypes)
                                WRITE(0,'(2a)') 'Line read in cfg file: ', Trim(phrase)
                                STOP '< ReadCfg >'
                        END IF
                        ! Try to read atom definition again
                        READ(inp,'(a500)',iostat=io) phrase
                        IF (io.EQ.0) THEN
                                IF (Present(nAux)) THEN
                                        READ(phrase,*,iostat=io) xc(1:3,i), aux(1:nAux,i)
                                ELSE
                                        READ(phrase,*,iostat=io) xc(1:3,i)
                                END IF
                        ELSE
                                WRITE(0,'(a,i0)') 'Error when reading in CFG file coordinates for atom' , i
                                WRITE(0,'(2a)') 'Line read in cfg file: ', Trim(phrase)
                                STOP '< ReadCfg >'
                        END IF
                END IF
        ELSE
                WRITE(0,'(a,i0)') 'Error when reading in CFG file coordinates for atom' , i
                WRITE(0,'(2a)') 'Line read in cfg file: ', Trim(phrase)
                STOP '< ReadCfg >'
        END IF
        iTyp(i) = nTypes
     END DO

     ! Real coordinates
     xp(:,1:im) = MatMul(at(:,:),xc(:,1:im))

     DEALLOCATE(xc)

  END SUBROUTINE ReadCfg

  SUBROUTINE WriteCfg(xp, iTyp, im, at, nTypes, mass, label, out_alat, out, &
        mask, nAux_int, aux_int, nAux_real, aux_real, aux_title)
    ! Write structure in output file connected to unit out using Cfg (atomeye) format
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written
    !  Additional columns are used to print auxiliary properties
    !   nAux_int, nAux_real: number of auxiliary properties to print (int and real types)
    !   aux_int(1:nAux_int), aux_real(1:nAux_real): corresponding property for atom i
    !   aux_title(1:nAux_int): title for auxiliary property with integer type
    !   aux_title(nAux_int+1:nAux_int+nAux_real): title for auxiliary property with real type

    USE mathStruct
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
    INTEGER, dimension(:), intent(in) :: iTyp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    INTEGER, intent(in) :: nTypes
    REAL(kind(0.d0)), dimension(:), intent(in) :: mass
    CHARACTER(len=5), dimension(:), intent(in) :: label
    REAL(kind(0.d0)), intent(in) :: out_alat
    INTEGER, intent(in) :: out
    LOGICAL, dimension(:), intent(in), optional :: mask
    INTEGER, intent(in), optional :: nAux_int, nAux_real
    INTEGER, dimension(:,:), optional :: aux_int
    REAL(kind(0.d0)), dimension(:,:), optional :: aux_real
    CHARACTER(len=50), dimension(:), intent(in), optional :: aux_title

    INTEGER :: i, ic, j, n
    CHARACTER(len=50) :: out_format
    LOGICAL, dimension(:), allocatable :: local_mask
    LOGICAL :: test_aux_int, test_aux_real, label_defined
    REAL(kind(0.d0)), dimension(:,:), allocatable :: xc
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
    REAL(kind(0.d0)) :: alat

    ! Check if atom types have been defined
    label_defined=.FALSE.
    DO n=1, nTypes
       IF (label(n).NE.aChar(n+64)) THEN
               label_defined=.true.
               EXIT
       END IF
    END DO
    IF (.NOT.label_defined) THEN
            WRITE(0,'(a)') 'WARNING < WriteCfg >'
            WRITE(0,'(a)') 'It seems that atom labels have not been defined'
            DO n=1, Max(nTypes,1)
               WRITE(0,'(a,i0,2a)') '  label for atom type ', n,': ', label(n)
            END DO
            WRITE(0,'(a)') 'You can define different labels using command label(1)="Fe" in Babel input file'
    END IF

    ! Check if lattice vectors have been defined and calculate reduced
    ! atom coordinates
    IF ( SUM( at(1:3,1:3)**2 ).LE.Distance_Zero2 ) THEN
            WRITE(0,'(a)') 'You need to define lattice vectors at(1:3,i) to be&
                & able to write configuration in cfg format'
            STOP '< WriteCfg >'
    END IF
    CALL Mat3Inv(at,inv_at)
    IF (Allocated(xc)) Deallocate(xc)
    ALLOCATE(xc(1:3,1:im))
    xc(1:3,1:im) = MatMul( inv_at(1:3,1:3), xp(1:3,1:im) )
    WHERE(abs(xc(1:3,1:im)).LE.1.d-15) xc(:,:)=0.d0     ! Problem with atomeye otherwhise with negative zero (-1.d-16)

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
                    STOP '< WriteCfg >'
            END IF
            IF (Size(aux_int,1).LT.nAux_int)  THEN
                    WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                    STOP '< WriteCfg >'
            END IF
            test_aux_int=.TRUE.
    ELSE
            test_aux_int=.FALSE.
    END IF
    IF ( Present(nAux_real) ) THEN
            IF ( (.NOT.Present(aux_real)) .OR. (.NOT.Present(aux_title)) ) THEN
                    WRITE(0,'(a)') 'All auxiliary properties have to be defined'
                    STOP '< WriteCfg >'
            END IF
            IF (Size(aux_real,1).LT.nAux_real)  THEN
                    WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                    STOP '< WriteCfg >'
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
    
    write(out,'(a,i0)')'Number of particles = ', Count(local_mask(1:im))
    write(out,'(a,g24.16,a)')'A = ', alat, ' Angstrom (basic length-scale)'
    do j=1,3
       write(out,'(a,i0)') '# Unit cell vector #', j
       do ic=1,3
          write(out,'(A,I1,A,I1,A,g24.16,A)')'H0(',j,',',ic,') = ',at(ic,j)/alat,' A'
       end do
    end do
    write(out,'(A)')'.NO_VELOCITY.'
    IF (test_aux_int.AND.test_aux_real) THEN
            write(out,'(a,i0)')'entry_count = ', 3+nAux_int+nAux_real
            DO i=1, nAux_int+nAux_real
               WRITE(out,'(a,i0,2a)') 'auxiliary[',i-1,'] = ', Trim(aux_title(i))
            END DO
            WRITE(out_format,'(a,3(i0,a))') &
                '(',3,'(g24.16,1x),',nAux_int,'(1x,i0),',nAux_real,'(1x,g24.16))'
            DO n=1, Max(nTypes,1)
               WRITE(out,'(f0.3)') mass(n)
               WRITE(out,'(a)') label(n)
               do i=1, im
                  IF (iTyp(i).NE.n) Cycle
                  IF (local_mask(i)) write(out,out_format) xc(1:3,i), aux_int(1:nAux_int,i), &
                           aux_real(1:nAux_real,i)
               end do
            END DO
    ELSEIF (test_aux_int) THEN
            write(out,'(a,i0)')'entry_count = ', 3+nAux_int
            DO i=1, nAux_int
               WRITE(out,'(a,i0,2a)') 'auxiliary[',i-1,'] = ', Trim(aux_title(i))
            END DO
            WRITE(out_format,'(a,2(i0,a))') &
                '(',3,'(g24.16,1x),',nAux_int,'(1x,i0)),'
            DO n=1, Max(nTypes,1)
               WRITE(out,'(f0.3)') mass(n)
               WRITE(out,'(a)') label(n)
               do i=1, im
                  IF (iTyp(i).NE.n) Cycle
                  IF (local_mask(i)) write(out,out_format) xc(1:3,i), aux_int(1:nAux_int,i)
               end do
            END DO
    ELSEIF (test_aux_real) THEN
            write(out,'(a,i0)')'entry_count = ', 3+nAux_real
            DO i=1, nAux_real
               WRITE(out,'(a,i0,2a)') 'auxiliary[',i-1,'] = ', Trim(aux_title(i))
            END DO
            WRITE(out_format,'(a,2(i0,a))') &
                '(',3,'(g24.16,1x),',nAux_real,'(1x,g24.16))'
            DO n=1, Max(nTypes,1)
               WRITE(out,'(f0.3)') mass(n)
               WRITE(out,'(a)') label(n)
               do i=1, im
                  IF (iTyp(i).NE.n) Cycle
                  IF (local_mask(i)) write(out,out_format) xc(1:3,i), aux_real(1:nAux_real,i)
               end do
            END DO
    ELSE
            write(out,'(a,i0)')'entry_count = ', 3
            WRITE(out_format,'(a,i0,a)') '(',3,'(g24.16,1x))'
            DO n=1, Max(nTypes,1)
               WRITE(out,'(f0.3)') mass(n)
               WRITE(out,'(a)') label(n)
               do i=1, im
                  IF (iTyp(i).NE.n) Cycle
                  IF (local_mask(i)) write(out,out_format) xc(1:3,i)
               end do           
            END DO
    ENDIF


    Deallocate(local_mask, xc)

    
  END SUBROUTINE WriteCfg
  
END MODULE cfg_module
