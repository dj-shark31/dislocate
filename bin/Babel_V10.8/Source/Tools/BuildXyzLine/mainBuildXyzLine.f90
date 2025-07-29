PROGRAM BuildXyzLine

  USE Structure_module
  IMPLICIT NONE

  INTEGER :: imm, im
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp
  INTEGER, dimension(:), allocatable :: iTyp
  REAL(kind(0.d0)), dimension(3,3) :: at
  REAL(kind(0.d0)) :: xMin, xMax, y0, z0, dx
  INTEGER :: i
  CHARACTER(len=100) :: outFile

  ! Input parameters
  WRITE(6,'(a)') "Number of atoms"
  READ(5,*) im
  WRITE(6,'(a)') "x min"
  READ(5,*) xMin
  WRITE(6,'(a)') "x max"
  READ(5,*) xMax
  WRITE(6,'(a)') "y0"
  READ(5,*) y0
  WRITE(6,'(a)') "z0"
  READ(5,*) z0
  WRITE(6,'(a)') "Name of the Xyz output file"
  READ(5,*) outFile

  ! Allocation and Initialization
  imm=im
  ALLOCATE(xp(1:3,1:imm))
  ALLOCATE(iTyp(1:imm))
  iTyp(1:im)=1          ! All atoms have type 1
  label(1)="A"          ! Default name for atom types
  at(1:3,1:3) = 0.d0

  ! Atom coordinates
  dx=(xMax-xMin)/dble(im-1)
  DO i=1, im
     xp(1,i) = xMin + dble(i-1)*dx
     xp(2,i) = y0
     xp(3,i) = z0
  END DO

  ! Write Xyz file
  CAll WriteStructure(xp, iTyp, im, at, outFile, "xyz")
  DEALLOCATE(xp, iTyp)

END PROGRAM BuildXyzLine
