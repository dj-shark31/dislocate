MODULE ndm_module
  ! Binary structure file used by NDM to backup atom positions (itesauv=)
  !  file *.coutposition

CONTAINS

  !===============================================

  FUNCTION GetImmNDM(inp) RESULT(imm)
    ! Get maximal number of atoms imm from *.coutposition file connected to unit inp
    
    IMPLICIT NONE
    INTEGER, intent(in) :: inp
    INTEGER :: imm

    INTEGER :: iCinType
    REAL(kind(0.d0)), dimension(3,3) :: at
    REAL(kind(0.d0)), dimension(3) :: zl

    imm=0
    READ(inp) iCinType
    IF (iCinType.GE.2) THEN
            READ(inp) at
    ELSE
            READ(inp) zl
    END IF
    READ(inp) imm
    REWIND(inp)

  END FUNCTION GetImmNDM

  !===============================================

  SUBROUTINE ReadNDM(xp, iTyp, im, at, nTypes, inp)

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(out) ::xp
    INTEGER, dimension(:), intent(out) :: iTyp
    INTEGER, intent(out) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: at
    INTEGER, intent(inout) :: nTypes
    INTEGER, intent(in) :: inp

    INTEGER :: iCinType, n
    REAL(kind(0.d0)), dimension(3) :: zl

    READ(inp) iCinType
    IF (iCinType.GE.2) THEN
            READ(inp) at
    ELSE
            READ(inp) zl
            at(:,:)=0.d0
            at(1,1)=zl(1) ; at(2,2)=zl(2) ; at(3,3)=zl(3)
    END IF
    READ(inp) im
    READ(inp) ityp(:)
    READ(inp) xp(:,:)

    ! Distances are in cm in NDM 
    ! I want them in angstroms

    at(:,:)=1.e8*at(:,:)
    xp(:,:)=1.e8*xp(:,:)

    ! Number of atom types
    n=MaxVal(ityp(1:im))
    nTypes=Max(nTypes,n)

  END SUBROUTINE ReadNDM

  !===============================================

  SUBROUTINE WriteNDM(xp, iTyp, im, at, out, mask)
    !  If mask is given, only atoms for which mask(i)=.TRUE. are written

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(in) ::xp
    INTEGER, dimension(:), intent(in) :: iTyp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    INTEGER, intent(in) :: out
    LOGICAL, dimension(:), intent(in), optional :: mask

    ! Local variables
    INTEGER :: i, im_new, imm
    REAL(kind(0.d0)), dimension(:,:), allocatable :: xp_new
    INTEGER, dimension(:), allocatable :: iTyp_new

    WRITE(out) 2        ! iCinType
    WRITE(out) at(:,:)

    IF (Present(mask)) THEN
            imm=Size(iTyp)
            ALLOCATE(iTyp_new(1:imm))
            ALLOCATE(xp_new(1:3,1:imm))
            im_new=0
            xp_new(:,:)=0.d0
            iTyp_new(:)=0
            DO i=1, im
               IF (mask(i)) THEN
                       im_new=im_new+1
                       iTyp_new(im_new) = iTyp(i)
                       xp_new(:,im_new) = xp(:,i)
               END IF
            END DO
            WRITE(out) im_new
            WRITE(out) iTyp_new(:)
            WRITE(out) xp_new(:,:)
            DEALLOCATE(iTyp_new)
            DEALLOCATE(xp_new)
    ELSE
            WRITE(out) im
            WRITE(out) iTyp(:)
            WRITE(out) xp(:,:)
    END IF

  END SUBROUTINE WriteNDM

END MODULE ndm_module
