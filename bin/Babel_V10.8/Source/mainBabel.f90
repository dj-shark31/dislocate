PROGRAM Babel

  USE babel_data
  USE math
  USE elasticity_ani
  USE Volterra_Module
  USE Disloc_elasticity_ani
  USE LineCouple_elasticity_ani  
  USE LineForce_elasticity_ani  
  USE DDipoleModule
  USE structure_module
  USE babel_lecture
  USE strain_module
  USE core_module
  USE slabModule
  USE loopModule
  USE neighboursModule
  USE Symmetry3DModule 
  IMPLICIT NONE
 
  ! Number of neighbours for each atom (input and reference structure)
  INTEGER, dimension(:), allocatable :: nNeigh, nNeigh0
  ! Neighbour indexes for each atom (input and reference structure)
  INTEGER, dimension(:,:), allocatable :: iNeigh, iNeigh0

  INTEGER, parameter :: verbosity_max=1

  INTEGER ::nAux_real, nAux_int, n, i, nAdd
  CHARACTER(len=50), dimension(:), allocatable :: aux_title
  REAL(kind(0.d0)), dimension(:,:), allocatable :: aux_real
  INTEGER, dimension(:,:), allocatable :: aux_int
  LOGICAL, dimension(:), allocatable :: coreAtom
  !INTEGER, EXTERNAL :: iArgc
  LOGICAL :: strained
  CHARACTER(len=9) :: fileType
  LOGICAL, dimension(:), allocatable :: keepCut
  REAL(kind(0.d0)), dimension(1:6) :: sVoigt
  REAL(kind(0.d0)), dimension(:,:), allocatable :: uTotal, uRandom, uStrain, uAdd
  REAL(kind(0.d0)), dimension(1:3) :: uMean


  program_name="babel"

  ! Read name of the input file
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF

  ! Read input parameters in file 'input.dat'
  CALL Read_Babel(6)
  IF (verbosity.GE.10) THEN
          WRITE(6,*)
          WRITE(6,'(a)') 'Input file succesfully read'
          WRITE(6,*)
  END IF


  ! Atom initial coordinates
  IF (Allocated(xp0)) DEALLOCATE(xp0)
  IF (Allocated(iTyp0)) DEALLOCATE(iTyp0)
  IF (Allocated(xp)) THEN
          ALLOCATE(xp0(1:3,1:imm))
          ALLOCATE(iTyp0(1:imm))
          xp0(:,1:im) = xp(:,1:im)
          iTyp0(1:im) = iTyp(1:im)
          im0=im
          at0(:,:) = at(:,:)
  END IF

  ! Atom total displacement
  IF (Allocated(uTotal)) Deallocate(uTotal)
  Allocate(uTotal(1:3,1:imm))
  uTotal(:,:)=0.d0

  ! ===================================================
  ! Duplicate unit cell
  IF (duplicate) THEN
          CALL duplicate_cell(lat, xp, iTyp, im, at, 6, verbose=verbosity.GE.3)
          CALL duplicate_cell(lat, xp0, iTyp0, im0, at0, 6, verbose=.FALSE.)
          duplicate=.FALSE.
          lat(1:3)=1
  END IF

  ! ===================================================
  ! Symmetrize atomic positions
  IF (symmetrize) THEN
          CALL Symmetrize3DPositions(xp, iTyp, im, at, symGroup, 6, verbose=verbosity.GE.3)
          symmetrize=.FALSE.
  END IF

  ! ===================================================
  ! Rotate unit cell
  IF (rotate) THEN

          ! ---- Atom positions and periodicity vectors ---
          CALL rotate_cell(rot, xp, im, at, 6, verbose=verbosity.GE.3) 
          IF (verbosity.GE.10) THEN
                  WRITE(6,*)
                  WRITE(6,'(a)') 'Input unit cell succesfully rotated'
                  WRITE(6,*)
          END IF
          CALL rotate_cell(rot, xp0, im, at0, 6, verbose=.FALSE.) 
          IF (verbosity.GE.10) THEN
                  WRITE(6,*)
                  WRITE(6,'(a)') 'Reference unit cell succesfully rotated'
                  WRITE(6,*)
          END IF

          ! ----- Line defects -----------------------
          IF (l_loop)     CALL RotateLoops(rot,6,verbose=verbosity.GE.3)
          IF (l_dislo)    CALL RotateDislo(rot,6,verbose=verbosity.GE.3)
          IF (l_DDipole)  CALL RotateDDipole(rot,6,verbose=verbosity.GE.3)
          IF (l_lineCouple) CALL RotateLineCouple(rot,6,verbose=verbosity.GE.3)
          IF (lineForce)  CALL RotateLineForce(rot,6,verbose=verbosity.GE.3)

          ! ---- Homogeneous strain ------------------
          IF (Strain) THEN 
                  eStrain   = MatMul( rot, MatMul( eStrain, Transpose(rot) ) )
                  IF (verbosity.GE.3) THEN
                          WRITE(6,'(a)') '  new homogeneous strain applied:'
                          DO i=1, 3
                              WRITE(6,'(a,3g14.6,a)') '      | ', eStrain(i,1:3), ' |'
                          END DO
                          WRITE(6,*)
                          IF (anisotropic_elasticity) THEN
                                  sVoigt = MatMul( CVoigt, (/ eStrain(1,1), eStrain(2,2), estrain(3,3), &
                                      eStrain(2,3)+eStrain(3,2), eStrain(1,3)+eStrain(3,1), eStrain(1,2)+eStrain(2,1) /) )
                                  WRITE(6,'(a)') '  Corresponding stress (GPa):'
                                  WRITE(6,'(a,3g14.6,a)') '      | ', sVoigt(1), sVoigt(6), sVoigt(5), ' |'
                                  WRITE(6,'(a,3g14.6,a)') '      | ', sVoigt(6), sVoigt(2), sVoigt(4), ' |'
                                  WRITE(6,'(a,3g14.6,a)') '      | ', sVoigt(5), sVoigt(4), sVoigt(3), ' |'
                                  WRITE(6,*)
                          END IF
                  END IF
          END IF 

          ! ---- Elastic dipole for impurity -------------
          IF (impurity) THEN 
                  pImpurity = MatMul( rot, MatMul( pImpurity, Transpose(rot) ) )
                  IF (verbosity.GE.3) THEN
                          WRITE(6,'(a)') '  new impurity elastic dipole:'
                          DO i=1, 3
                              WRITE(6,'(a,3g14.6,a)') '      | ', pImpurity(i,1:3), ' |'
                          END DO
                          WRITE(6,*)
                  END IF
          END IF 

          !! ---- Translation vectors ----------
          !IF (translate) THEN
                  !uTranslate = MatMul( rot, uTranslate )
                  !IF (verbosity.GE.3) WRITE(6,'(a,3g14.6)') '  new translation vector: ', uTranslate(1:3)
          !END IF

          ! ---- Elastic constants ------------
          IF (anisotropic_elasticity) THEN
                  CVoigt = Rotate_CVoigt(CVoigt,rot)
                  CALL MatInv(CVoigt, inv_CVoigt)
                  IF (verbosity.GE.3) THEN
                          WRITE(6,'(a)') '  Rotated elastic constants (GPa)'
                          CALL Print_CVoigt(CVoigt, 6)
                          WRITE(6,'(a)') '  Inverse matrix (GPa^-1)'
                          CALL Print_CVoigt(inv_CVoigt,6)
                          WRITE(6,*)
                  END IF
          END IF

          rotate=.FALSE.
  END IF

  ! ===================================================
  ! Translate unit cell
  ! ----------------
  IF (translate) THEN
          CALL Translate_cell(uTranslate,xp,im,6,verbose=verbosity.GE.3)
          CALL Translate_cell(uTranslate,xp0,im0,6,verbose=verbosity.GE.3)
          IF (l_loop) CALL TranslateLoops(uTranslate,6,verbose=verbosity.GE.3)
          IF (l_dislo) CALL TranslateDislo(uTranslate,6,verbose=verbosity.GE.3)
          IF (l_DDipole) CALL TranslateDDipole(uTranslate,6,verbose=verbosity.GE.3)
          IF (lineForce)  CALL TranslateLineForce(uTranslate,6,verbose=verbosity.GE.3)
          IF (l_lineCouple) CALL TranslateLineCouple(uTranslate,6,verbose=verbosity.GE.3)
          translate=.false.
  END IF

  ! Extract a slab from unit cell
  IF (slab) THEN
          CALL MakeSlab(xp, iTyp, im, at, 6, verbosity.GE.3)
          slab=.FALSE.
  END IF
            
  ! Remove Volterra cut created by dislocation or loops
  IF (remove_cut) THEN
          nAdd = 0
          IF ((l_dislo.OR.l_loop.OR.l_DDipole)) THEN
                  ALLOCATE(keepCut(1:imm))
                  ALLOCATE(uAdd(1:3,1:imm))
                  IF (l_dislo) THEN

                          ! Version with both vacancy and interstitial loops
                          CALL RemoveAddDisloCut(xp0, xp, iTyp, im, keepCut, uAdd, at, 6, verbose=verbosity.GE.4) 
                          nAdd = nAdd + im-im0
                          iTyp0(:)=iTyp
                          im0=im
                          keep(:)=keep(:).AND.keepCut(:)
                          uTotal(:,:) = uTotal(:,:) + uAdd(:,:)
                  END IF
                  IF (l_DDipole) THEN

                          ! Version with both vacancy and interstitial loops
                          CALL RemoveAddDDipoleCut(xp0, xp, iTyp, im, keepCut, uAdd, at, 6, verbose=verbosity.GE.4)
                          nAdd = nAdd + im-im0
                          iTyp0(:)=iTyp
                          im0=im
                          keep(:)=keep(:).AND.keepCut(:)
                          uTotal(:,:) = uTotal(:,:) + uAdd(:,:)
                  END IF
                  IF (l_loop) THEN

                          ! Version with both vacancy and interstitial loops
                          CALL RemoveAddLoopCut(xp0, xp, iTyp, im, keepCut, uAdd, at, 6, verbose=verbosity.GE.4)
                          nAdd = nAdd + im-im0
                          iTyp0(:)=iTyp
                          im0=im
                          keep(:)=keep(:).AND.keepCut(:)
                          uTotal(:,:) = uTotal(:,:) + uAdd(:,:)
                  END IF
                  DEALLOCATE(keepCut)
                  DEALLOCATE(uAdd)
          END IF
  END IF

  strained=(l_loop.OR.l_dislo.OR.l_DDipole.OR.LineForce.OR.l_lineCouple.OR.strain.OR.read_uFile)

  ! Create dislocations, line-forces and line-force couples
  IF (strained) THEN
          IF (Allocated(uStrain)) Deallocate(uStrain)
          Allocate(uStrain(1:3,1:imm))
          CALL elastic_strain(uStrain, 6)
          uTotal(:,:) = uTotal(:,:) + uStrain(:,:)
          DEALLOCATE(uStrain)
          l_loop=.false.
          l_dislo=.false.
          l_DDipole=.false.
          lineForce=.false.
          l_lineCouple=.false.
          strain=.false.
          read_uFile=.false.
  END IF

  ! Add noise to atomic positions keeping fixed the gravity center
  IF (xNoise.GT.0.d0) THEN
          ALLOCATE(uRandom(3,im))
          CALL Random_Number( uRandom )
          uRandom(1:3,1:im) = 2.d0*xNoise*( uRandom(1:3,1:im) - 0.5d0)
          uMean(1:3) = Sum( uRandom(1:3,1:im), 2)/dble(im)
          DO i=1, im
             uRandom(1:3,i) = uRandom(1:3,i) - uMean(1:3)
          END DO
          uTotal(1:3,1:im) = uTotal(1:3,1:im) + uRandom(1:3,1:im)
          DEALLOCATE(uRandom)
  END IF


  ! Add displacement
  IF (.NOT.initial) THEN
          ! initial=.FALSE. => apply displacement to atom coordinates
          xp(:,1:im) = xp0(:,1:im) + uTotal(:,1:im)              ! New cartesian coordinates
  ELSE
          xp(:,1:im) = xp0(:,1:im)
          at=at0
  END IF

  ! Clip atoms in unit cell applying periodic boundary conditions
  IF (clipAtom) THEN
          CALL Clip_Atoms(xp, im, at, 6, verbosity.GE.3)
          CALL Clip_Atoms(xp0, im0, at0, 6, verbosity.GE.3)
          clipAtom=.false.
  END IF

  nAux_real=0 ; nAux_int=0; n=0
  IF (out_neighbours) nAux_int=nAux_int+1
  IF (out_id)          nAux_int=nAux_int+1
  IF (out_core)        nAux_int=nAux_int+1
  IF (out_x) nAux_real=nAux_real+1
  IF (out_y) nAux_real=nAux_real+1
  IF (out_z) nAux_real=nAux_real+1
  IF (strained.AND.out_displacement)  nAux_real=nAux_real+3
  IF (strained.AND.out_pressure)      nAux_real=nAux_real+1
  IF (strained.AND.out_VonMises)      nAux_real=nAux_real+1
  IF (strained.AND.out_stress)        nAux_real=nAux_real+6
  IF (strained.AND.out_elasticStrain) nAux_real=nAux_real+6
  IF (strained.AND.out_Ebinding)      nAux_real=nAux_real+1
  IF ( (nAux_int.NE.0) .OR. (nAux_real.NE.0) ) THEN
          ALLOCATE(aux_title(1:nAux_real+nAux_int))
          IF (nAux_int.NE.0) ALLOCATE(aux_int(1:nAux_int,1:imm))
          IF (nAux_real.NE.0) ALLOCATE(aux_real(1:nAux_real,1:imm))
          IF (out_neighbours) THEN      ! Print number of neighbours for each atom
                  CALL InitNeighbours(imm)
                  ALLOCATE(nNeigh(1:imm))               ; nNeigh(:) = 0
                  ALLOCATE(iNeigh(1:max_nNeigh, 1:imm)) ; iNeigh(:,:) = 0
                  CALL BuildNeighbours(im, xp, at, nNeigh, iNeigh, 6)
                  aux_int(n+1,:) = nNeigh(:)
                  IF (outLammps) THEN
                          aux_title(n+1)='nNeigh'
                  ELSE
                          aux_title(n+1)='Number of neighbours'
                  END IF
                  n=n+1
          END IF
          IF (out_id) THEN      ! Print atom id
                  DO i=1, imm
                     aux_int(n+1,i) = i
                  END DO
                  IF (outLammps) THEN
                          aux_title(n+1)='id'
                  ELSE
                          aux_title(n+1)='Atom id'
                  END IF
                  n=n+1
          END IF
          IF (out_core) THEN ! Print atom core index
                  ALLOCATE(CoreAtom(1:imm))
                  CALL Core(CoreAtom, 6) 
                  aux_int(n+1,:)=0
                  WHERE (CoreAtom(:)) aux_int(n+1,:)=1
                  IF (outLammps) THEN
                          aux_title(n+1)='core'
                  ELSE
                          aux_title(n+1)='Core index'
                  END IF
                  n=n+1
                  DEALLOCATE(CoreAtom)
                
          END IF
          IF (out_x) THEN       ! x coordinate
                  aux_real(n-nAux_int+1,:) = xp(1,:)
                  IF (outLammps) THEN
                          aux_title(n+1)='x'
                  ELSE
                          aux_title(n+1)='Atom x coordinate (A)'
                  END IF
                  n=n+1
          END IF
          IF (out_y) THEN       ! x coordinate
                  aux_real(n-nAux_int+1,:) = xp(2,:)
                  IF (outLammps) THEN
                          aux_title(n+1)='y'
                  ELSE
                          aux_title(n+1)='Atom y coordinate (A)'
                  END IF
                  n=n+1
          END IF
          IF (out_z) THEN       ! x coordinate
                  aux_real(n-nAux_int+1,:) = xp(3,:)
                  IF (outLammps) THEN
                          aux_title(n+1)='z'
                  ELSE
                          aux_title(n+1)='Atom z coordinate (A)'
                  END IF
                  n=n+1
          END IF
          IF (strained.AND.out_displacement) THEN  ! Print atom displacement
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = uTotal(1:3,:)
                  IF (outLammps) THEN
                          aux_title(n+1)='ux'
                          aux_title(n+2)='uy'
                          aux_title(n+3)='uz'
                  ELSE
                          aux_title(n+1)='Total displacement Ux'
                          aux_title(n+2)='                   Uy'
                          aux_title(n+3)='                   Uz (A)'
                  END IF
                  n=n+3
          END IF
          IF (strained.AND.out_pressure) THEN      ! Print pressure
                  aux_real(n-nAux_int+1,:) = pressure(:)
                  IF (outLammps) THEN
                          aux_title(n+1)='p'
                  ELSE
                          aux_title(n+1) = 'Pressure (GPa)'
                  END IF
                  n=n+1
          END IF
          IF (strained.AND.out_VonMises) THEN      ! Print Von-Misès equivalent shear stress
                  aux_real(n-nAux_int+1,:) = VonMises(:)
                  IF (outLammps) THEN
                          aux_title(n+1)='sVonMises'
                  ELSE
                          aux_title(n+1) = 'Von-Misès equivalent shear stress (GPa)'
                  END IF
                  n=n+1
          END IF
          IF (strained.AND.out_stress) THEN        ! Print atomic stress (6 components)
                  aux_real(n-nAux_int+1:n-nAux_int+6,:)=stress(1:6,:)
                  IF (outLammps) THEN
                          aux_title(n+1)='sxx'
                          aux_title(n+2)='syy'
                          aux_title(n+3)='szz'
                          aux_title(n+4)='syz'
                          aux_title(n+5)='sxz'
                          aux_title(n+6)='sxy'
                  ELSE
                          aux_title(n+1) = 'Stress Sxx'
                          aux_title(n+2) = '       Syy'
                          aux_title(n+3) = '       Szz'
                          aux_title(n+4) = '       Syz'
                          aux_title(n+5) = '       Sxz'
                          aux_title(n+6) = '       Sxy (GPa)'
                  END IF
                  n=n+6
          END IF
          IF (strained.AND.out_elasticStrain) THEN        ! Print atomic strain (6 components)
                  aux_real(n-nAux_int+1:n-nAux_int+6,:)=elasticStrain(1:6,:)
                  IF (outLammps) THEN
                          aux_title(n+1)='exx'
                          aux_title(n+2)='eyy'
                          aux_title(n+3)='ezz'
                          aux_title(n+4)='eyz'
                          aux_title(n+5)='exz'
                          aux_title(n+6)='exy'
                  ELSE
                          aux_title(n+1) = 'Strain Exx'
                          aux_title(n+2) = '       Eyy'
                          aux_title(n+3) = '       Ezz'
                          aux_title(n+4) = '     2*Eyz'
                          aux_title(n+5) = '     2*Exz'
                          aux_title(n+6) = '     2*Exy'
                  END IF
                  n=n+6
          END IF
          IF (strained.AND.out_Ebinding) THEN     ! Print interaction energy with impurity
                  aux_real(n-nAux_int+1,:) = Ebinding(:)
                  IF (outLammps) THEN
                          aux_title(n+1)='eBind'
                  ELSE
                          aux_title(n+1) = 'Ebinding (eV)'
                  END IF
                  n=n+1
          END IF
  END IF

  IF (n.NE.nAux_real+nAux_int) THEN
          WRITE(0,'(a,i0)') ' n         = ', n
          WRITE(0,'(a,i0)') ' nAux_int = ', nAux_int
          WRITE(0,'(a,i0)') ' nAux_real = ', nAux_real
  END IF

  ! ==== Output structure ========================
  IF (verbosity.GE.verbosity_max) THEN
          WRITE(6,*)
          IF (outXyz)       WRITE(6,'(2a)') 'Write output structure with Xyz format in file ', outFile
          IF (outOnlyAtoms) WRITE(6,'(2a)') 'Write one line per atom in file ', outFile
          IF (outCfg)       WRITE(6,'(2a)') 'Write output structure with Cfg format in file ', outFile
          IF (outGin)       WRITE(6,'(2a)') 'Write output structure with Gin format in file ', outFile
          IF (outLisa)      WRITE(6,'(2a)') 'Write output structure with Lisa format in file ', outFile
          IF (outSiesta)    WRITE(6,'(2a)') 'Write output structure with Siesta format in file ', outFile
          IF (outNDM)       WRITE(6,'(2a)') 'Write output structure with binary NDM format in file ', outFile
          IF (outLammps)    WRITE(6,'(2a)') 'Write output structure with Lammps dump format in file ', outFile
          IF (outPoscar)    WRITE(6,'(2a)') 'Write output structure with Poscar dump format in file ', outFile
          IF (strained.AND.initial) THEN
                  WRITE(6,'(a)') '  Initial atom coordinates are kept'
                  IF (at_defined) &
                        WRITE(6,'(a)') '  Initial periodicity vectors are kept'
          ELSE IF (strained) THEN
                  WRITE(6,'(a)') '  Displacement field &
                      &is added to initial atom coordinates'
                  IF (at_defined) &
                        WRITE(6,'(a)') '  Equivalent homogeneous strain is&
                                & applied to periodicity vectors'
          END IF
          IF (remove_cut) THEN
                  n=Count(Keep(1:im))
                  WRITE(6,'(2x,3(i0,a))') im-n, ' atoms have been removed and ', & 
                           nAdd, ' atoms have been inserted during Volterra process :  ', &
                        n, ' atoms now'
          END IF
          IF (outXyz) THEN
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': atom label'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': x coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': y coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 4, ': z coordinate (A)'
                  DO n=1, nAux_int+nAux_real
                     WRITE(6,'(a,i2,2a)') '   column ', n+4, ': ', aux_title(n)
                  END DO
          END IF
          IF (outOnlyAtoms) THEN
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': x coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': y coordinate (A)'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': z coordinate (A)'
                  DO n=1, nAux_int+nAux_real
                     WRITE(6,'(a,i2,2a)') '   column ', n+3, ': ', aux_title(n)
                  END DO
          END IF
          IF (outCfg) THEN
                  DO n=1, nAux_int+nAux_real
                     WRITE(6,'(a,i2,2a)') '   property ', n-1, ': ', aux_title(n)
                  END DO
          END IF
          IF (outLammps) THEN
                  WRITE(6,'(a,i2,a)') '   column ', 1, ': atom id'
                  WRITE(6,'(a,i2,a)') '   column ', 2, ': atom type'
                  WRITE(6,'(a,i2,a)') '   column ', 3, ': reduced coordinate xsu'
                  WRITE(6,'(a,i2,a)') '   column ', 4, ':                    ysu'
                  WRITE(6,'(a,i2,a)') '   column ', 5, ':                    zsu'
                  DO n=1, nAux_int+nAux_real
                     WRITE(6,'(a,i2,2a)') '   column ', n+5, ': ', aux_title(n)
                  END DO
          END IF
          WRITE(6,*)
  END IF

  IF (outXyz.OR.outOnlyAtoms.OR.outCfg.OR.outGin.OR.outLisa.OR.outSiesta.OR.outNDM.OR.outLammps.OR.outPoscar) THEN
          IF (outXyz) fileType='xyz'
          IF (outOnlyAtoms) fileType='onlyAtoms'
          IF (outCfg) fileType='cfg'
          IF (outGin) fileType='gin'
          IF (outSiesta) fileType='siesta'
          IF (outLisa) fileType='lisa'
          IF (outNDM) fileType='ndm'
          IF (outLammps) fileType='lammps'
          IF (outPoscar) fileType='poscar'
          IF ( (nAux_int.GT.0).AND.(nAux_real.GT.0) ) THEN
                  CALL WriteStructure(xp, iTyp, im, at, outFile, fileType, mask=keep(:), &
                        nAux_int=nAux_int, aux_int=aux_int(:,:), &
                        nAux_real=nAux_real, aux_real=aux_real(:,:), &
                        aux_title=aux_title(:))
          ELSE IF (nAux_int.GT.0) THEN
                  CALL WriteStructure(xp, iTyp, im, at, outFile, fileType, mask=keep(:), &
                        nAux_int=nAux_int, aux_int=aux_int(:,:), &
                        aux_title=aux_title(:))
          ELSE IF (nAux_real.GT.0) THEN
                  CALL WriteStructure(xp, iTyp, im, at, outFile, fileType, mask=keep(:), &
                        nAux_real=nAux_real, aux_real=aux_real(:,:), &
                        aux_title=aux_title(:))
          ELSE
                  CALL WriteStructure(xp, iTyp, im, at, outFile, fileType, mask=keep(:))
          END IF
  END IF

  ! ======== Reference structure =================================

  IF (refXyz.OR.refCfg.OR.refGin.OR.refLisa.OR.refSiesta.OR.refNDM) THEN

          IF (verbosity.GE.verbosity_max) THEN
                  WRITE(6,*)
                  IF (refXyz) WRITE(6,'(2a)')    'Write reference structure with Xyz format in file ', refFile
                  IF (refCfg) WRITE(6,'(2a)')    'Write reference structure with Cfg format in file ', refFile
                  IF (refGin) WRITE(6,'(2a)')    'Write reference structure with Gin format in file ', refFile
                  IF (refLisa) WRITE(6,'(2a)')   'Write reference structure with Lisa format in file ', refFile
                  IF (refSiesta) WRITE(6,'(2a)') 'Write reference structure with Siesta format in file ', refFile
                  IF (refNDM) WRITE(6,'(2a)')    'Write reference structure with binary NDM format in file ', refFile
                  IF (refLammps) WRITE(6,'(2a)') 'Write reference structure with Lammps format in file ', refFile
                  IF (refPoscar) WRITE(6,'(2a)') 'Write reference structure with Poscar format in file ', refFile
                  WRITE(6,'(a)') '  Initial atom coordinates are kept'
                  IF (at_defined) &
                        WRITE(6,'(a)') '  Initial periodicity vectors are kept'
                  IF (remove_cut) THEN
                          n=Count(Keep(1:im))
                          WRITE(6,'(2(2x,i0,a))') im-n, ' atoms have been removed during Volterra process:', &
                                n, ' atoms left'
                  END IF
                  IF (refXyz) THEN
                          WRITE(6,'(a,i2,a)') '   column ', 1, ': atom label'
                          WRITE(6,'(a,i2,a)') '   column ', 2, ': x coordinate (A)'
                          WRITE(6,'(a,i2,a)') '   column ', 3, ': y coordinate (A)'
                          WRITE(6,'(a,i2,a)') '   column ', 4, ': z coordinate (A)'
                  END IF
                  WRITE(6,*)
          END IF

          IF (refXyz)    fileType='xyz'
          IF (refCfg)    fileType='cfg'
          IF (refGin)    fileType='gin'
          IF (refSiesta) fileType='siesta'
          IF (refLisa)   fileType='lisa'
          IF (refNDM)    fileType='ndm'
          IF (refLammps) fileType='lammps'
          IF (refLammps) fileType='poscar'
          CALL WriteStructure(xp0, iTyp, im, at0, refFile, fileType, mask=keep(:))


  END IF

END PROGRAM Babel
