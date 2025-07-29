MODULE associate_data

  IMPLICIT NONE
  SAVE

  INTEGER :: verbosity

  REAL(kind(0.d0)), parameter :: distance_zero=1.d-4
  REAL(kind(0.d0)), parameter :: distance_zero2=distance_zero**2

  !========================================================
  ! Periodicity vectors for input and reference structure
  REAL(kind(0.d0)), dimension(1:3, 1:3) :: at, at0
  ! Maximal number of atoms in simulation box (used to allocate tables)
  INTEGER :: imm=0
  ! Real number of atoms in simulation box (input and reference structures)
  INTEGER :: im=0, im0=0
  ! Atom real coordinates: xp(1:3,:) (input and reference structures)
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp, xp0
  ! Atom type (input and reference structures)
  INTEGER, dimension(:), allocatable :: iTyp, iTyp0

  ! Atom index in input and reference structures
  !  atom i in input structure corresponds to atom ind(i) in reference structure
  !  atom i in reference structure corresponds to atom ind0(i) in input structure
  INTEGER, dimension(:), allocatable :: ind, ind0


  !========================================================
  ! Format of the structure files (output and reference)
  LOGICAL :: outXyz, outCfg, outGin, outLisa, outSiesta
  ! Output and input structure files
  CHARACTER(len=100) :: outFile
  ! Sort atoms according to their index in output structure file
  LOGICAL :: out_sort
  ! Print atom index in output structure file
  LOGICAL :: out_index

  !========================================================
  ! Apply or not periodic boundary conditions when looking for closest replica
  LOGICAL :: periodic

END MODULE associate_data
