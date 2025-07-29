MODULE Lisa_module
  ! Structure file in format used by Lisa in her Siesta calculations

CONTAINS

  !=============================================================

  FUNCTION GetImmLisa(inp) RESULT(imm)
    ! Get maximal number of atoms imm from lisa file connected to unit inp
    
    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    INTEGER :: imm

    REAL(kind(0.d0)), dimension(1:3) :: R
    CHARACTER(len=5) :: symbol
    INTEGER :: i, io

    ! First pass to determine number of atoms
    CALL Comment(inp)
    READ(inp,*) ! Lattice parameter
    DO i=1, 3
       CALL Comment(inp)
       READ(inp,*) ! Periodicity vectors
    END DO
    imm=0
    DO 
       CALL Comment(inp)
       READ(inp,*,iostat=io) R(1:3), i, symbol
       IF (io.NE.0) Exit
       imm=imm+1
    END DO
    Rewind(inp)

  END FUNCTION GetImmLisa

  !=============================================================

  SUBROUTINE ReadLisa(xp, iTyp, im, at, nTypes, label, inp)
    ! Read structure in input file connected to inp

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(out) ::xp
    INTEGER, dimension(:), intent(out) :: iTyp
    INTEGER, intent(out) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: at
    INTEGER, intent(inout) :: nTypes
    CHARACTER(len=5), dimension(:), intent(out) :: label
    INTEGER, intent(in) :: inp

    REAL(kind(0.d0)) :: alat
    REAL(kind(0.d0)), dimension(1:3) :: R
    CHARACTER(len=5) :: symbol
    INTEGER :: i, io, n

    ! Initialization
    alat = 1.d0
    xp(:,:) = 0.d0
    iTyp(:) = 0
    im = 0
    at(:,:) = 0.d0

    ! Read lattice parameter
    CALL Comment(inp)
    READ(inp,*) alat ! Lattice parameter is defined

    ! Read periodicity vectors
    DO i=1, 3
       CALL Comment(inp)
       READ(inp,*) at(1:3,i)
    END DO

    ! Read atom definitions
    io=0
    DO 
       CALL Comment(inp)
       READ(inp,*,iostat=io) R(1:3), i, symbol
       IF (io.NE.0) Exit
       im=im+1
       IF ( (im.GT.size(xp,2)).OR.(im.GT.size(iTyp,1)) ) THEN
               WRITE(0,'(a)') 'Too many atoms found in input structure file'
               WRITE(0,'(a)') 'You need to increase imm in input file'
               STOP '< ReadLisa >' 
       END IF
       xp(:,im) = R(:)
       iTyp(im)=i
       IF (i.GT.Size(label(:),1)) THEN
                  WRITE(0,'(a)') '  Too many atom types found => atom labels may be wrong'
                  WRITE(0,'(a)') '  You need to increase parameter max_nTypes in structure.f90 and to recompile'
                  WRITE(0,'(a,i0)') '    current value: ', Size(label(:),1)
                  STOP '< ReadLisa >'
       END IF
       label(i)=symbol
    END DO

    n=MaxVal(ityp(1:im))
    nTypes=Max(nTypes,n)

    ! Distance normalization
    xp(:,1:im) = alat*xp(:,1:im)
    at(:,:) = alat*at(:,:)

  END SUBROUTINE ReadLisa

  !=============================================================

  SUBROUTINE WriteLisa(xp, iTyp, im, at, label, out, mask)
    ! Write structure in output file connected to unit out using Lisa format
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
    INTEGER, dimension(:), intent(in) :: iTyp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    CHARACTER(len=5), dimension(:), intent(in) :: label
    INTEGER, intent(in) :: out
    LOGICAL, dimension(:), intent(in), optional :: mask

    ! Local variables
    REAL(kind(0.d0)), parameter :: alat=1.d0
    INTEGER :: i

    ! Lattice parameter
    WRITE(out,'(g24.16)') alat
    WRITE(out,*)

    ! Periodicity vectors
    DO i=1, 3
       WRITE(out,'(3g24.16)') at(1:3,i) ! Basis vector i
    END DO
    WRITE(out,*)

    ! Atom
    IF (Present(mask)) THEN
            DO i=1, im
               IF (mask(i)) WRITE(out,'(3(g24.16,1x),i0,1x,a,1x,i0)') &
                        xp(1:3,i), iTyp(i), label(iTyp(i)), i
            END DO
    ELSE
            DO i=1, im
               WRITE(out,'(3(g24.16,1x),i0,1x,a,1x,i0)') xp(1:3,i), iTyp(i), label(iTyp(i)), i
            END DO
    END IF

  END SUBROUTINE WriteLisa

END MODULE Lisa_module
