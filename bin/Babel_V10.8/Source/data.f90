MODULE babel_data

  USE Symmetry3DModule, ONLY : sym3D_group_t
  SAVE

  CHARACTER(len=20) :: program_name     ! babel, inter

  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter :: Distance_Zero2=Distance_Zero*Distance_Zero

  ! Identity matrix
  REAL(kind(0.d0)), dimension(1:3,1:3), parameter :: matId=ReShape( &
        (/ 1.d0, 0.d0, 0.d0, &
           0.d0, 1.d0, 0.d0, &
           0.d0, 0.d0, 1.d0 /), (/ 3,3 /) )

  INTEGER :: verbosity
  INTEGER, parameter :: verbosity_debug=10
  LOGICAL :: debug

  ! Cutoff radius for elastic calculations
  REAL(kind(0.d0)) :: rc

  LOGICAL :: remove_cut ! If .true., remove atoms when creating dislocation dipoles

  !========================================================
  ! Keep fixed or not gravity center when creating line-defects
  LOGICAL :: fixGravity

  !========================================================
  ! Homogeneous strain applied and corresponding stress
  LOGICAL :: Strain
  REAL(kind(0.d0)), dimension(1:3,1:3) :: eStrain
  REAL(kind(0.d0)), dimension(1:6) :: sStrainVoigt
  LOGICAL :: induced_homogeneous_strain ! If true, the homogeneous strain induced by line-defects is considered

  !========================================================
  ! Solid displacement
  LOGICAL :: Translate
  REAL(kind(0.d0)), dimension(1:3) :: uTranslate

  !========================================================
  ! Strain due to impurity: elastic dipole (in eV)
  LOGICAL :: impurity
  REAL(kind(0.d0)), dimension(1:3,1:3) :: pImpurity

  !========================================================
  REAL(kind(0.d0)) :: factorE   ! Factor to apply to energy

  !========================================================
  ! Periodic images
  !  if xImages is true, consider nxImages images in direction
  !       at(1:3,1) and -at(1:3,1)
  LOGICAL :: xImages, yImages, zImages
  INTEGER :: nxImages, nyImages, nzImages
  LOGICAL :: clipDisplacement   ! If .true., apply periodic boundary conditions to atom displacements
  LOGICAL :: clipAtom   ! If .true., apply periodic boundary conditions to atom coordinates

  !========================================================
  ! Control of properties printed on output
  LOGICAL :: out_neighbours     ! Print number of neighbours for each atom
  LOGICAL :: out_id             ! Print atom id
  LOGICAL :: out_pattern        ! Print pattern index for each atom
  LOGICAL :: out_x, out_y, out_z! Print atom x, y or z coordinates
  LOGICAL :: out_pressure       ! Print pressure if true
  LOGICAL :: out_VonMises       ! Print equivalent Von-Misès shear stress if true
  LOGICAL :: out_displacement   ! Print atom displacement if true
  LOGICAL :: out_stress         ! Print atomic stress (6 components)
  LOGICAL :: out_strain         ! Print atomic total strain (6 components)
  LOGICAL :: out_elasticStrain         ! Print atomic elastic strain (6 components)
  LOGICAL :: out_plasticStrain         ! Print atomic plastic strain (6 components)
  LOGICAL :: out_gradDisplacement ! Print gradient of total displacement (9 components)
  LOGICAL :: out_gradElasticDisplacement ! Print gradient of elastic displacement (9 components)
  LOGICAL :: out_Nye            ! Print Nye tensor (9 components)
  LOGICAL :: out_BurgersDensity ! Print density of dislocation defined by bNye(1:3) and lNye(1:3)
  REAL(kind(0.d0)), dimension(1:3) :: bNye, lNye ! Burgers and line direction for dislocation density
  LOGICAL :: out_Ebinding ! Print interaction energy with impurity strain
  LOGICAL :: out_core     ! Print index 1 (0) for atoms belonging (not) to line-defect cores

  !========================================================
  ! Elastic constants in Voigt notation (GPa) and inverse matrix
  LOGICAL :: anisotropic_elasticity, isotropic_elasticity
  REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt, inv_CVoigt
  REAL(kind(0.d0)) :: CVoigt_noise      ! noise to add for creating line-defects

  !========================================================
  ! If duplicate=.true. the unit cell is duplicated lat(1), lat(2)... times in each
  ! direction given by at(1:3,1), at(1:3,2)...
  LOGICAL :: duplicate
  INTEGER, dimension(3) :: lat

  !========================================================
  ! If rotate=.true. the unit cell is rotated
  !   rot(:,:) is the rotation matrix
  LOGICAL :: rotate
  REAL(kind(0.d0)), dimension(3,3) :: rot

  !========================================================
  ! If symmetrize=.true. the structure is symmetrized
  !   with symmetry operations defined by sym
  LOGICAL :: symmetrize
  TYPE(sym3D_group_t), SAVE :: symGroup


  !========================================================
  ! Format of the structure files (output and reference)
  LOGICAL :: outXyz, outCfg, outGin, outLisa, outSiesta, outNDM, outOnlyAtoms, outLammps, outPoscar
  LOGICAL :: refXyz, refCfg, refGin, refLisa, refSiesta, refNDM, refLammps, refPoscar
  ! Output and input structure files
  CHARACTER(len=100) :: outFile, refFile, inpFile

  !========================================================
  ! If initial=.true., keep initial coordinates when creating a dislocation
  !   else add dislocation displacement field to atom coordinates
  LOGICAL :: initial

  !========================================================
  ! DEFINITION OF THE STRUCTURE
  ! Lattice vector coordinates (A): at(1:3,1), at(2:3,2), ...
  ! (input and reference structures)
  LOGICAL :: at_defined, at0_defined
  REAL(kind(0.d0)), dimension(1:3, 1:3) :: at, at0, at1
  ! Maximal number of atoms in simulation box (used to allocate tables)
  INTEGER :: imm=0
  ! Real number of atoms in simulation box (input and reference structures)
  INTEGER :: im=0, im0=0, im1=0
  ! Atom real coordinates: xp(1:3,:) (input and reference structures)
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp, xp0, xp1
  ! Atom type (input and reference structures)
  INTEGER, dimension(:), allocatable :: iTyp, iTyp0, iTyp1
  ! keep(i): if .false., atom i is removed
  LOGICAL, dimension(:), allocatable :: keep

  !========================================================
  ! Displacement read/write in a file 
  ! uFile(1:3,i) displacement of atom i
  REAL(kind(0.d0)), dimension(:,:), allocatable :: uFile
  LOGICAL :: read_uFile

  !========================================================
  ! Control of Eulerian coordinates
  ! x = x0 + u(x), where u(x) is the displacement calculated in x
  INTEGER :: max_Euler  ! Maximal number of cycles in Eulerian scheme
  ! max_Euler=0 => Lagrangian coordinates
  REAL(kind(0.d0)) :: delta_Euler       ! Absolute tolerance for displacement convergency

  !========================================================
  ! Noise to add to atomic positions
  REAL(kind(0.d0)) :: xNoise

END MODULE babel_data
