SUBROUTINE Version(out)

  USE babel_data
  IMPLICIT NONE

  INTEGER, intent(in) :: out
  CHARACTER(len=8) :: date
  CHARACTER(len=10) :: time
  INTEGER, parameter :: verbosity_max=1

  IF (verbosity.LT.verbosity_max) RETURN

  CALL Date_and_Time(date, time)

  WRITE(out,*)
  WRITE(out,'(3a)')  '  Program ', Program_Name, ' (V10.8 - June 2023)'
  WRITE(out,'(a)')   '        (c) SRMP - CEA Saclay - France'
  WRITE(out,'(a)')   '   contact: emmanuel.clouet@cea.fr'
  WRITE(out,*)
  WRITE(out,'(12a)') '  Execution time: ', &
        date(1:4), '-', date(5:6), '-', date(7:8), '  ', &
        time(1:2), ':', time(3:4), ':', time(5:10)
  WRITE(out,*)

END SUBROUTINE Version

