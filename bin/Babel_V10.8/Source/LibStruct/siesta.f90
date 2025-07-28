MODULE Siesta_module
  ! Structure file in format used by Siesta
  ! This corresponds to the format of the file *.STRUCT_IN and *.STRUCT_OUT read
  ! and written by Siesta

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter, private :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter, private :: Distance_Zero2=Distance_Zero*Distance_Zero

CONTAINS

  !=============================================================

  FUNCTION GetImmSiesta(inp) RESULT(imm)
    ! Get maximal number of atoms imm from siesta file connected to unit inp
    
    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    INTEGER :: imm

    INTEGER :: i

    imm=0

    ! Read periodicity vectors
    DO i=1, 3
       CALL Comment(inp)
       READ(inp,*) 
    END DO

    ! Read number of atoms
    CALL Comment(inp)
    READ(inp,*) imm

    REWIND(inp)

  END FUNCTION GetImmSiesta

  !=============================================================

  SUBROUTINE ReadSiesta(xp, iTyp, im, at, nTypes, inp)
    ! Read structure in input file connected to inp
    !   xp(1:3,i): Carthesian coordinates for atom i
    !   iTyp(i): type for atom i
    !   im: number of atoms
    !   at(1:3,i): Carthesian coordinates for periodicity vector i
    !   nTypes: number of different atom types

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(out) ::xp
    INTEGER, dimension(:), intent(out) :: iTyp
    INTEGER, intent(out) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: at
    INTEGER, intent(inout) :: nTypes
    INTEGER, intent(in) :: inp

    REAL(kind(0.d0)), dimension(1:3) :: xc
    INTEGER :: i, iDummy, n
    
    ! Initialization
    im = 0
    at(:,:) = 0.d0

    ! Read periodicity vectors
    DO i=1, 3
       CALL Comment(inp)
       READ(inp,*) at(1:3,i)
    END DO

    ! Read number of atoms
    CALL Comment(inp)
    READ(inp,*) im

    ! Maximal number of atoms in simulation box
    IF ( (im.GT.size(xp,2)).OR.(im.GT.size(iTyp,1)) ) THEN
             WRITE(0,'(a,i0)') '  dimension of xp(1:3,:): ', size(xp,2)
             WRITE(0,'(a,i0)') '  dimension of iTyp(:): ', size(iTyp,1)
             WRITE(0,'(a,i0)') '  number of atoms read in Siesta file: ', im
             STOP '< ReadSiesta >'
    END IF
    
    ! Initialization
    xp(:,:) = 0.d0
    iTyp(:) = 0

    ! Read atom definitions
    DO i=1, im
       CALL Comment(inp)
       READ(inp,*) iTyp(i), iDummy, xc(1:3)
       xp(1:3,i) = MatMul( at(:,:), xc(:) )
    END DO

    ! Number of atom types
    n=MaxVal(ityp(1:im))
    nTypes=Max(nTypes,n)

  END SUBROUTINE ReadSiesta

  !=============================================================

  SUBROUTINE WriteSiesta(xp, iTyp, im, at, out, mask)
    ! Write structure in output file connected to unit out using Siesta format
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written

    USE MathStruct   
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
    INTEGER, dimension(:), intent(in) :: iTyp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    INTEGER, intent(in) :: out
    LOGICAL, dimension(:), intent(in), optional :: mask

    ! Local variables
    REAL(kind(0.d0)), dimension(1:3) :: xc
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
    INTEGER :: i

    ! Check that periodicity vectors are defined
    IF ( SUM( at(1:3,1:3)**2 ).LE.Distance_Zero2 ) THEN
            WRITE(0,'(a)') 'You need to define lattice vectors at(1:3,i) to be&
                & able to use Siesta format'
            STOP '< WriteSiesta >'
    END IF

    ! Inverse matrix of periodicity vectors
    CALL Mat3Inv(at,inv_at)

    ! Periodicity vectors
    DO i=1, 3
       WRITE(out,'(3g24.16)') at(1:3,i) ! Basis vector i
    END DO

    ! Number of atoms
    IF (Present(mask)) THEN
            WRITE(out,'(i0)') Count(mask(1:im))
    ELSE
            WRITE(out,'(i0)') im
    END IF

    ! Atom
    IF (Present(mask)) THEN
            DO i=1, im
               xc(:) = MatMul( inv_at(:,:), xp(:,i) )
               IF (mask(i)) WRITE(out,'(2(i0,1x),3(g24.16))') &
                        iTyp(i), 0, xc(1:3)
            END DO
    ELSE
            DO i=1, im
               xc(:) = MatMul( inv_at(:,:), xp(:,i) )
               WRITE(out,'(2(i0,1x),3(g24.16))') &
                        iTyp(i), 0, xc(1:3)
            END DO
    END IF

  END SUBROUTINE WriteSiesta

END MODULE Siesta_module
