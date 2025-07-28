MODULE StringModule

CONTAINS

  FUNCTION stringUpperCase(strIn) RESULT(strOut)
    ! Convert string strIn in upper case

    IMPLICIT NONE
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j
    INTEGER, parameter :: jMin=iachar("a"), jMax=iaChar("z")

    do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>=jmin .and. j<=jMax ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
    end do

  END FUNCTION stringUpperCase

  FUNCTION stringLowerCase(strIn) RESULT(strOut)
    ! Convert string strIn in lower case

    IMPLICIT NONE
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j
    INTEGER, parameter :: jMin=iachar("A"), jMax=iaChar("Z")

    do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>=jmin .and. j<=jMax ) then
               strOut(i:i) = achar(iachar(strIn(i:i))+32)
          else
               strOut(i:i) = strIn(i:i)
          end if
    end do

  END FUNCTION stringLowerCase

END MODULE StringModule
