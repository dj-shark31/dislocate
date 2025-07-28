SUBROUTINE WriteInput(inpFile)

  IMPLICIT NONE
  CHARACTER(len=*), intent(in) :: inpFile

  CHARACTER(len=200) :: line
  INTEGER :: io

  OPEN(file=inpFile,unit=50,action='write')
  DO
     READ(5,'(a)',iostat=io) line
     IF (io.NE.0) EXIT
     WRITE(50,'(a)') Trim(line)
  END DO
  CLOSE(50)

END SUBROUTINE WriteInput
