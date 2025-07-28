PROGRAM PrepareDrag

  USE Math
  USE babel_data
  USE structure_module
  USE drag_lecture
  IMPLICIT NONE
 
  INTEGER, parameter :: verbosity_max=1
  REAL(kind(0.d0)), parameter :: Zero=1.d-6
  REAL(kind(0.d0)), parameter :: Zero2=Zero*Zero

  INTEGER ::i, cUnit
  REAL(kind(0.d0)), dimension(:,:), allocatable :: uTotal, ds, xc, xc0, xc1
  REAL(kind(0.d0)), dimension(1:3,1:3) :: at_temp, at0_temp, at1_temp, inv_at0, inv_at1
  !INTEGER, EXTERNAL :: iArgc
  CHARACTER(len=6) :: fileType
  REAL(kind(0.d0)), dimension(:,:), allocatable :: uRandom, dg
  REAL(kind(0.d0)), dimension(1:3) :: uMean
  REAL(kind(0.d0)) :: inv_dg, lambda

  program_name="prepareDrag"

  ! Read name of the input file
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF

  ! Read input parameters in file 'input.dat'
  CALL Read_Drag(6)

  ALLOCATE(uTotal(1:3,1:imm))   ! Total displacement (cartesian)
  ALLOCATE(ds(1:3,1:imm))       ! Total displacement (reduced)
  Allocate(xp(1:3,1:imm))       ! Cartesian coordinates in interpolated structure
  Allocate(xc(1:3,1:imm))       ! Reduced coordinates in interpolated structure
  ALLOCATE(xc1(1:3,1:imm))      ! Reduced coordinates in final structure
  ALLOCATE(xc0(1:3,1:imm))      ! Reduced coordinates in initial structure

  IF ( (.NOT.clipDisplacement).AND.(.NOT.reduced) ) THEN
          ! without periodic boundary conditions
          uTotal(:,1:im) = xp1(:,1:im) - xp0(:,1:im)
          at(:,:) = at0(:,:) + zeta*( at1(:,:) - at0(:,:) )
          xp(:,1:im) = xp0(:,1:im) + zeta*uTotal(:,1:im)
  ELSE
          ! with periodic boundary conditions
          IF (.NOT.at_defined) THEN
                  WRITE(0,'(a)') 'Basis vectors at(1:3,1:3) undefined'
                  WRITE(0,'(a)') '  => could not apply periodic boundary conditions to atom displacement'
                  WRITE(0,'(a)') '     You need to define at(1:3,1:3) or to set clipDisplacement and reduced to .false. '
                  STOP '< mainPrepareDrag >'
          END IF

          IF ( (MatNorm2( at1(:,:) - at0(:,:) ).GT.Zero).AND.(reduced.EQV..FALSE.) ) THEN
                  WRITE(0,'(a)') 'Basis vectors in initial and final configurations differ'
                  WRITE(0,'(a)') '  => use "reduced=.true." to constraint reduced coordinates'
                  WRITE(0,'(a)') '     instead of cartesian coordinates'
                  STOP '< mainPrepareDrag >'
          END IF

          ! Modify periodicity vectors to take into account the fact that unit cells can
          ! have been duplicated in each direction
          at0_temp(1:3,1) = at0(1:3,1)/dble(nxSlab)
          at0_temp(1:3,2) = at0(1:3,2)/dble(nySlab)
          at0_temp(1:3,3) = at0(1:3,3)/dble(nzSlab)
          CALL MatInv(at0_temp,inv_at0)

          at1_temp(1:3,1) = at1(1:3,1)/dble(nxSlab)
          at1_temp(1:3,2) = at1(1:3,2)/dble(nySlab)
          at1_temp(1:3,3) = at1(1:3,3)/dble(nzSlab)
          CALL MatInv(at1_temp,inv_at1)

          ! Reduced coordinate of the final state
          xc1(:,1:im)=MatMul(inv_at1(:,:),xp1(:,1:im))
          
          ! Reduced coordinates of the initial state
          xc0(:,1:im)=MatMul(inv_at0(:,:),xp0(:,1:im))

          ! Displacement between initial and final state
          ds(:,1:im) = xc1(:,1:im) - xc0(:,1:im)
          IF (clipDisplacement) ds(:,1:im) = ds(:,1:im) - aNInt(ds(:,1:im))
          at_temp(:,:) = at0_temp(:,:) + zeta*( at1_temp(:,:) - at0_temp(:,:) )
          uTotal(:,1:im) = MatMul( at_temp(:,:), ds(:,1:im) )

          ! Interpolate structure 
          IF (reduced) THEN
                  xc(:,1:im) = xc0(:,1:im) + zeta*ds(:,1:im)
                  xp(:,1:im) = MatMul( at_temp(:,:), xc(:,1:im) )
          ELSE
                  xp(:,1:im) = xp0(:,1:im) + zeta*uTotal(:,1:im)
          END IF

          ! Interpolate periodicity vectors
          at(:,:) = at0(:,:) + zeta*( at1(:,:) - at0(:,:) )

          DEALLOCATE(xc0) ; DEALLOCATE(xc1)

  END IF

  ! Write displacement in constraint file
  IF (constrFile.NE.'') THEN
          IF (constrFile.EQ.'-') THEN
                  cUnit=6
          ELSE
                  cUnit=62
                  OPEN(file=constrFile, unit=cUnit, action='write')
          END IF
          IF (reduced) THEN
                  IF (matNorm2(ds).LE.Zero*matNorm2(0.5*(xc0+xc1))) THEN
                          WRITE(0,'(a)') 'Displacement between initial and final states is null'
                          STOP '< PrepareDrag >'
                  END IF
                  DO i=1, im
                     WRITE(cUnit,'(3(g24.16,1x))') ds(1:3, i)
                  END DO 
                  WRITE(cUnit,'(a)') 'REDUCED :  constraint is applied on reduced coordinates'
          ELSE
                  IF (matNorm2(uTotal).LE.Zero*matNorm2(0.5*(xp0+xp1))) THEN
                          WRITE(0,'(a)') 'Displacement between initial and final states is null'
                          STOP '< PrepareDrag >'
                  END IF
                  DO i=1, im
                     WRITE(cUnit,'(3(g24.16,1x))') uTotal(1:3, i)
                  END DO
                  WRITE(cUnit,'(a)') 'CARTESIAN :  constraint is applied on cartesian coordinates'
          END IF
          IF (cUnit.NE.6) CLOSE(cUnit)
  END IF


  ! Add noise to atomic positions keeping fixed the gravity center
  IF (xNoise.GT.0.d0) THEN
          ! Generate random displacement
          ALLOCATE(uRandom(3,im)) ; ALLOCATE(dg(3,im))
          CALL Random_Number( uRandom )
          uRandom(1:3,1:im) = 2.d0*xNoise*( uRandom(1:3,1:im) - 0.5d0)
          ! Remove solid displacement
          uMean(1:3) = Sum( uRandom(1:3,1:im), 2)/dble(im)
          DO i=1, im
             uRandom(1:3,i) = uRandom(1:3,i) - uMean(1:3)
          END DO
          ! Project displacement in hyperplane
          inv_dg = 1.d0 / Sqrt( Sum( uTotal(1:3,1:im)**2 ) )
          dg(1:3,1:im) = uTotal(1:3,1:im)*inv_dg
          lambda = Sum( uRandom(1:3,1:im)*dg(1:3,1:im) )
          uRandom(1:3,1:im) = uRandom(1:3,1:im) - lambda*dg(1:3,1:im)
          ! Apply random displacement
          xp(1:3,1:im) = xp(1:3,1:im) + uRandom(1:3,1:im)
          DEALLOCATE(uRandom) ; DEALLOCATE(dg)
  END IF


  ! Clip atoms in unit cell applying periodic boundary conditions
  IF (clipAtom) THEN
          CALL Clip_Atoms(xp, im, at, 6, verbosity.GE.3)
          clipAtom=.false.
  END IF

  IF (verbosity.GE.verbosity_max) THEN
          WRITE(6,*)
          IF (outXyz) WRITE(6,'(2a)') 'Write output structure with Xyz format in file ', outFile
          IF (outCfg) WRITE(6,'(2a)') 'Write output structure with Cfg format in file ', outFile
          IF (outLammps)    WRITE(6,'(2a)') 'Write output structure with Lammps dump format in file ', outFile
          IF (outPoscar)    WRITE(6,'(2a)') 'Write output structure with Poscar dump format in file ', outFile
          IF (outXyz) THEN
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': atom label'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': x coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': y coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 4, ': z coordinate (A)'
          END IF
          IF (outCfg) THEN
                  WRITE(6,'(a,i2,a)')  '   property ', 0, ': atom type'
          END IF
          IF (outLammps) THEN
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': atom id'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': atom type'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': reduced coordinate xsu'
                  WRITE(6,'(a,i2,a)') '   column ', 4, ':                    ysu'
                  WRITE(6,'(a,i2,a)') '   column ', 5, ':                    zsu'
          END IF
          WRITE(6,*)
  END IF

  !======================================================================
  ! Output structure
  IF (outXyz.OR.outCfg.OR.outGin.OR.outLisa.OR.outSiesta.OR.outNDM.OR.outLammps.OR.outPoscar) THEN
          IF (outXyz) fileType='xyz'
          IF (outCfg) fileType='cfg'
          IF (outGin) fileType='gin'
          IF (outSiesta) fileType='siesta'
          IF (outLisa) fileType='lisa'
          IF (outNDM) fileType='ndm'
          IF (outLammps) fileType='lammps'
          IF (outPoscar) fileType='poscar'
          CALL WriteStructure(xp, iTyp, im, at, outFile, fileType)
  END IF

END PROGRAM PrepareDrag
