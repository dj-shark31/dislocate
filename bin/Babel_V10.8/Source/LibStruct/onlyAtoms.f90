MODULE onlyAtoms_module
  ! only atoms are printed with one line per atom and associated property

CONTAINS

  !=============================================================

  SUBROUTINE WriteOnlyAtoms(xp, im, out_alat, out, &
        mask, nAux_int, aux_int, nAux_real, aux_real, aux_title)
    ! Write structure in output file connected to unit out using only atom coordinates
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written
    !  Additional columns are used to print auxiliary properties
    !   nAux_int, nAux_real: number of auxiliary properties to print (int and real types)
    !   aux_int(1:nAux_int), aux_real(1:nAux_real): corresponding property for atom i
    !   aux_title(1:nAux_int): title for auxiliary property with integer type
    !   aux_title(nAux_int+1:nAux_int+nAux_real): title for auxiliary property with real type

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), intent(in) :: out_alat
    INTEGER, intent(in) :: out
    LOGICAL, dimension(:), intent(in), optional :: mask
    INTEGER, intent(in), optional :: nAux_int, nAux_real
    INTEGER, dimension(:,:), optional :: aux_int
    REAL(kind(0.d0)), dimension(:,:), optional :: aux_real
    CHARACTER(len=50), dimension(:), intent(in), optional :: aux_title

    INTEGER :: i
    REAL(kind(0.d0)) :: alat
    CHARACTER(len=50) :: out_format
    LOGICAL, dimension(:), allocatable :: local_mask
    LOGICAL :: test_aux_int, test_aux_real

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
                    STOP '< WriteOnlyAtoms >'
            END IF
            IF (Size(aux_int,1).LT.nAux_int)  THEN
                    WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                    STOP '< WriteOnlyAtoms >'
            END IF
            test_aux_int=.TRUE.
    ELSE
            test_aux_int=.FALSE.
    END IF
    IF ( Present(nAux_real) ) THEN
            IF ( (.NOT.Present(aux_real)) .OR. (.NOT.Present(aux_title)) ) THEN
                    WRITE(0,'(a)') 'All auxiliary properties have to be defined'
                    STOP '< WriteOnlyAtoms >'
            END IF
            IF (Size(aux_real,1).LT.nAux_real)  THEN
                    WRITE(0,'(a)') 'Problem with size of auxiliary properties'
                    STOP '< WriteOnlyAtoms >'
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


    IF (test_aux_int.AND.test_aux_real) THEN
            WRITE(out_format,'(a,i0,a)') '(a,', nAux_int+nAux_real,'(", ",a))'
            WRITE(out,out_format) '# x, y, z', &
                (Trim(AdjustL(aux_title(i))),i=1,nAux_int+nAux_real) ! Title
            WRITE(out_format,'(a,3(i0,a))') &
                '(',3,'(g24.16,1x),',nAux_int,'(1x,i0),',nAux_real,'(1x,g24.16))'
            DO i=1, im
               IF (local_mask(i)) WRITE(out,out_format) xp(1:3,i)/alat, &
                        aux_int(1:nAux_int,i),  aux_real(1:nAux_real,i)
            END DO
    ELSEIF (test_aux_int) THEN
            WRITE(out_format,'(a,i0,a)') '(a,', nAux_int,'(", ",a))'
            WRITE(out,out_format) '# x, y, z', &
                (Trim(AdjustL(aux_title(i))),i=1,nAux_int) ! Title
            WRITE(out_format,'(a,2(i0,a))') &
                '(',3,'(g24.16,1x),',nAux_int,'(1x,i0))'
            DO i=1, im
               IF (local_mask(i)) WRITE(out,out_format) xp(1:3,i)/alat, &
                        aux_int(1:nAux_int,i)
            END DO
    ELSEIF (test_aux_real) THEN
            WRITE(out_format,'(a,i0,a)') '(a,', nAux_real,'(", ",a))'
            WRITE(out,out_format) '# x, y, z', &
                (Trim(AdjustL(aux_title(i))),i=1,nAux_real) ! Title
            WRITE(out_format,'(a,2(i0,a))') &
                '(',3,'(g24.16,1x),',nAux_real,'(1x,g24.16))'
            DO i=1, im
               IF (local_mask(i)) WRITE(out,out_format) xp(1:3,i)/alat, &
                        aux_real(1:nAux_real,i)
            END DO
    ELSE
            WRITE(out,'(a)') '# x, y, z'                ! Title
            WRITE(out_format,'(a,i0,a)') '(',3,'(g24.16,1x))'
            DO i=1, im
               IF (local_mask(i)) WRITE(out,out_format) xp(1:3,i)/alat
            END DO
    END IF

  END SUBROUTINE WriteOnlyAtoms

END MODULE onlyAtoms_module
