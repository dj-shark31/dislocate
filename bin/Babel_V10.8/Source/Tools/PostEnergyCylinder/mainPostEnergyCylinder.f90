PROGRAM PostEnergyCylinder

  USE structure_module
  USE cfg_module
  USE mathStruct
  IMPLICIT NONE

  ! Atom real coordinates: xp(1:3,:) (input and reference structures)
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xp
  ! Lattice vector coordinates (A): at(1:3,1), at(2:3,2), ...
  REAL(kind(0.d0)), dimension(1:3, 1:3) :: at
  ! Maximal number of atoms in simulation box (used to allocate tables)
  INTEGER :: imm
  ! Real number of atoms in simulation box
  INTEGER :: im
  ! Atom type (input and reference structures)
  INTEGER, dimension(:), allocatable :: iTyp
  REAL(kind(0.d0)), dimension(:,:), allocatable :: aux
  !!$CHARACTER(len=50), dimension(:), allocatable :: aux_title
  REAL(kind(0.d0)), dimension(:), allocatable :: eAtom, dAtom
  REAL(kind(0.d0)), dimension(1:3) :: cCylinder, lCylinder
  REAL(kind(0.d0)), dimension(1:3) :: dx, ds
  REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
  REAL(kind(0.d0)) :: r0, r1, dr, Ecylinder, Eref
  INTEGER :: i, nAtom, dnAtom
  CHARACTER(len=100) :: input_file, outFile

  ! Input file
  WRITE(6,'(a)') 'Name of the cfg file:'
  READ(5,*) input_file
  OPEN(file=input_file, unit=50, action='read', status='old')
  imm=GetImmCfg(50)
  ALLOCATE(xp(1:3,1:imm))
  ALLOCATE(iTyp(1:imm))
  ALLOCATE(aux(1,1:imm))
  CALL ReadCfg(xp, iTyp, im, at, nTypes, mass, label, 50, nAux=1, aux=aux)
  CLOSE(50)

  ! Table containing energy per atom
  Allocate(Eatom(1:imm))
  Eatom(1:im) = aux(1,1:im)
  Eatom(im+1:imm) = 0.d0
  DEALLOCATE(aux)

  ! Inverse matrix of lattice vectors
  CALL Mat3Inv(at, inv_at)

  ! Definition of the cylinder
  WRITE(6,'(a)') 'Cylinder axis origin' 
  READ(5,*) cCylinder(1:3)
  WRITE(6,'(a)') 'Cylinder axis direction'
  READ(5,*) lCylinder(1:3)
  lCylinder(1:3) = lCylinder(1:3)/Sqrt( Sum(lCylinder(1:3)**2) )

  ! Atom distance to cylinder axis
  Allocate(dAtom(1:imm))
  DO i=1, im
    dx(:) = xp(:,i) - cCylinder(:)
    ! Apply periodic boundary conditions
      ds(:) = MatMul(inv_at(:,:), dx(:) )
      ds(:) = ds(:) - aNInt( ds(:) )
      dx = MatMul( at(:,:), ds(:) )
    dAtom(i) = Sqrt( &
         ( dx(2) * lCylinder(3) - dx(3) * lCylinder(2) )**2 &
       + ( dx(3) * lCylinder(1) - dx(1) * lCylinder(3) )**2 &
       + ( dx(1) * lCylinder(2) - dx(2) * lCylinder(1) )**2 )
  END DO

  !=== DEBUG =====
  !!$ALLOCATE(aux_title(2))
  !!$Allocate(aux(2,1:imm))
  !!$aux(1,1:im)=Eatom(1:im)
  !!$aux(2,1:im)=dAtom(1:im)
  !!$aux_title(1)='Energy per atom'
  !!$aux_title(2)='distance'
  !!$OPEN(file='temp.cfg', unit=61, action='write')
  !!$CALL WriteCfg(xp, iTyp, im, at, 61,  &
        !!$nAux_real=2, aux_real=aux(:,:), aux_title=aux_title(:))
  !!$Deallocate(aux_title)
  !!$Deallocate(aux)
  !!$CLOSE(61)
  !!$STOP
  !=== DEBUG =====

  ! Radius increment
  WRITE(6,'(a)') 'Radius increment'
  READ(5,*) dr

  ! Reference energy
  WRITE(6,'(a)') 'Energy per atom in perfect crystal'
  READ(5,*) Eref

  ! Output file
  WRITE(6,'(a)') 'Output file'
  READ(5,*) outFile
  OPEN(file=outFile, unit=60, action='write', status='unknown')
  WRITE(60,'(a,g20.6)') '# radius increment: dr = ', dr
  WRITE(60,'(a,g30.16)') '# energy per atom in perfect crystal: Eref = ', Eref
  WRITE(60,'(a,3(g20.6,1x))') '# cylinder axis origin: ', cCylinder(1:3)
  WRITE(60,'(a,3(g20.6,1x))') '# cylinder axis direction: ', lCylinder(1:3)
  WRITE(60,'(2a)') '# name of the input file: ', input_file
  WRITE(60,'(a,i0)'), '# number of atoms: ', im
  WRITE(60,*)
  WRITE(60,'(a)') '# r, nAtom, E(r), dE(r)'

  ! Integrate energy
  nAtom=0
  Ecylinder=0.d0
  r1=0.d0
  DO 
     r0=r1 ; r1=r0+dr
     dnAtom=0
     !!$WRITE(6,*) ! DEBUG
     !!$WRITE(6,'(a,g20.6)') 'r = ', r0+0.5*dr     ! DEBUG
     DO i=1, im
        IF ( (dAtom(i).GT.r0).AND.(dAtom(i).LE.r1) ) THEN
                ECylinder = ECylinder + eAtom(i)
                dnAtom = dnAtom+1
                !!$WRITE(6,'(a,i0)') 'atom ', i    ! DEBUG
        END IF
     END DO
     !!$READ(5,*)  ! DEBUG
     IF (dnAtom.NE.0) THEN
             nAtom = nAtom + dnAtom
             WRITE(60,'(g20.6,1x,i0,2(g20.6,1x))') &
                r0+0.5d0*dr, nAtom, ECylinder, ECylinder-dble(nAtom)*Eref
     END IF
     IF (nAtom.GE.im) Exit
  END DO

  CLOSE(60)
  
  Deallocate(Eatom)

END PROGRAM PostEnergyCylinder
