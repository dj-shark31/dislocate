PROGRAM Displacement

  USE babel_data
  USE structure_module
  USE displacementModule
  USE neighboursModule
  USE gradDisplacementModule
  USE gradElasticDisplacementModule
  USE strainFromDisplacementModule
  USE nyeTensorModule 

  IMPLICIT NONE
 
  INTEGER, parameter :: verbosity_max=1

  ! Total displacement for each atom
  REAL(kind(0.d0)), dimension(:,:), allocatable :: uTotal

  ! Total and elastic strain for each atom (Voigt notation)
  REAL(kind(0.d0)), dimension(:,:), allocatable :: eAtomicStrain, eAtomicElasticStrain, eAtomicPlasticStrain

  ! Number of neighbours for each atom (input and reference structure)
  INTEGER, dimension(:), allocatable :: nNeigh, nNeigh0
  ! Neighbour indexes for each atom (input and reference structure)
  INTEGER, dimension(:,:), allocatable :: iNeigh, iNeigh0

  INTEGER :: nAux_real, nAux_int, n, i, iX, iY, iZ
  CHARACTER(len=50), dimension(:), allocatable :: aux_title
  REAL(kind(0.d0)), dimension(:,:), allocatable :: aux_real
  INTEGER, dimension(:,:), allocatable :: aux_int  
  REAL(kind(0.d0)), dimension(1:3) :: uSolid, u2Solid 
  REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at0, inv_at, epsi
  REAL(kind(0.d0)), dimension(:,:), allocatable :: xc, xc0
  !INTEGER, EXTERNAL :: iArgc
  REAL(kind(0.d0)) :: u2, uTemp2, bNorm, lNorm, dU2
  REAL(kind(0.d0)), dimension(1:3) :: uNew, uTemp
  CHARACTER(len=9) :: fileType
  LOGICAL :: ok

  program_name="displacement"

  ! Read name of the input file
  IF (iArgc().LE.0) THEN
          WRITE(6,'(a)') 'Name of the input file'
          READ(5,*) input_file
  ELSE
          CALL getArg(1,input_file)
  END IF

  ! Read input parameters in file 'input.dat'
  CALL Read_displacement(6)

  IF (out_BurgersDensity) THEN
          lNorm = Sqrt( Sum( lNye(1:3)**2 ) )
          bNorm = Sqrt( Sum( bNye(1:3)**2 ) )
          IF ( ( lNorm.LE.distance_zero ).OR.( bNorm.LE.distance_zero ) ) THEN
                  WRITE(0,'(a)') 'You need to define a direction for the line &
                          and the Burgers vector for projecting dislocation density'
                  WRITE(0,'(a,3g14.6)') ' lNye(1:3) = ', lNye(1:3)
                  WRITE(0,'(a,3g14.6)') ' bNye(1:3) = ', bNye(1:3)
                  STOP '< mainDisplacement >'
          END IF
  END IF

  ! Allocation
  ALLOCATE(uTotal(1:3,1:imm))

  ! Rotation
  IF (rotate) THEN  
          ! Rotate atom cartesian coordinates
          IF (verbosity.GE.verbosity_max) THEN
                  WRITE(6,'(a)') ' Rotation applied'
                  WRITE(6,'(a,3g14.6,a)') '   | ', rot(1,1:3), ' |'
                  WRITE(6,'(a,3g14.6,a)') '   | ', rot(2,1:3), ' |'
                  WRITE(6,'(a,3g14.6,a)') '   | ', rot(3,1:3), ' |'
                  WRITE(6,*)
          END IF
          CALL rotate_cell(rot, xp, im, at, 6, verbose=verbosity.GE.3) 
          CALL rotate_cell(rot, xp0, im0, at0, 6, verbose=verbosity.GE.3) 
          rotate=.FALSE.
  END IF

  ! Duplicate unit cell
  IF (duplicate) THEN
          IF (at_defined)  CALL duplicate_cell(lat, xp, iTyp, im, at, 6, verbose=verbosity.GE.3)
          IF (at0_defined) CALL duplicate_cell(lat, xp0, iTyp0, im0, at0, 6, verbose=.FALSE.)
          duplicate=.FALSE.
          Lat(1:3)=1
  END IF

  ! Total displacement
  IF (im0.LE.0) THEN
          ! Reference structure is not defined
          uTotal(:,:) = 0.d0
          IF (out_displacement) THEN
                  WRITE(6,*)
                  WRITE(6,'(a)') 'WARNING !!!!!!'
                  WRITE(6,'(a)') '  You need to define a reference structure to calculate total displacement'
                  WRITE(6,'(a)') '  ==> out_displacement set to .false.'    
                  out_displacement=.false.
          END IF
          IF (out_strain.OR.out_plasticStrain.OR.out_gradDisplacement) THEN
                  WRITE(6,*)
                  WRITE(6,'(a)') 'WARNING !!!!!!'
                  WRITE(6,'(a)') '  You need to define a reference structure to calculate total or plastic strains'
                  WRITE(6,'(a)') '     or gradient of total displacement'
                  WRITE(6,'(a)') '  ==> out_strain, out_plasticStrain and out_gradDisplacement set to .false.'    
                  out_strain=.false.
                  out_plasticStrain=.false.
                  out_gradDisplacement=.false.
          END IF
  ELSE
          IF (.NOT.clipDisplacement) THEN
                  ! without periodic boundary conditions
                  uTotal(:,1:im) = xp(:,1:im) - xp0(:,1:im)
          ELSE
                  ! with periodic boundary conditions
                  IF (.NOT.at_defined) THEN
                          WRITE(0,'(a)') 'Basis vectors at(1:3,1:3) undefined'
                          WRITE(0,'(a)') '  => could not apply periodic boundary conditions to atom displacement'
                          WRITE(0,'(a)') '     You need to define at(1:3,1:3) or to set clipDisplacement to .false. '
                          STOP '< mainDisplacement >'
                  END IF
                  CALL Mat3Inv(at,inv_at)
                  IF (Allocated(xc)) Deallocate(xc)
                  ALLOCATE(xc(1:3,1:imm))
                  xc(:,1:im)=MatMul(inv_at(:,:),xp(:,1:im))
                  IF (Allocated(xc0)) Deallocate(xc0)
                  CALL Mat3Inv(at0,inv_at0)
                  ALLOCATE(xc0(1:3,1:imm))
                  xc0(:,1:im)=MatMul(inv_at0(:,:),xp0(:,1:im))
                  xc(:,1:im) = xc(:,1:im) - aNInt(xc(:,1:im)-xc0(:,1:im))
                  uTotal(:,1:im) = MatMul(at(:,:),xc(:,1:im)) - MatMul(at0(:,:),xc0(:,1:im))

                  ! This loop is necessary when periodicity vectors are not orthogonal
                  DO i=1, im
                     u2 = Sum( uTotal(:,i)**2 )
                     uNew(:) = uTotal(:,i)
                     DO iX=-1,1; DO iY=-1,1 ; DO iZ=-1,1
                        IF ( (iX.EQ.0).AND.(iY.EQ.0).AND.(iZ.EQ.0) ) Cycle
                        uTemp(:) = uTotal(:,i) + dble(iX)*at(1:3,1) + dble(iY)*at(1:3,2) + dble(iZ)*at(1:3,3)
                        uTemp2 = Sum( uTemp(:)**2 )
                        IF ( uTemp2.LT.u2 ) THEN
                                u2 = uTemp2
                                uNew(:) = uTemp(:)
                        END IF
                     END DO; END DO; END DO
                     uTotal(:,i) = uNew(:)
                  END DO
                  DEALLOCATE(xc,xc0)
          END IF
  END IF

  ! Calculate solid displacement and remove it to keep fixed the gravity center
  IF (fixGravity) THEN
          uSolid(1:3) = Sum( uTotal(1:3,1:im), 2)/dble(im)
          DO i=1, im
            uTotal(1:3,i) = uTotal(1:3,i) - uSolid(1:3)
            xp(1:3,i) = xp(1:3,i) - uSolid(1:3)
          END DO
          IF (verbosity.GE.verbosity_max) THEN
                  WRITE(6,'(a,3g14.6)') 'Solid displacement removed &
                      &to keep fixed the gravity center: u(1:3) = ', uSolid(1:3)
                  WRITE(6,*)
          END IF
  END IF

  ! Add solid displacement
  IF (translate) THEN
          DO i=1, im
            uTotal(1:3,i) = uTotal(1:3,i) + uTranslate(1:3)
            xp(1:3,i) = xp(1:3,i) + uTranslate(1:3)
            xp0(1:3,i) = xp0(1:3,i) + uTranslate(1:3)
          END DO
          IF (verbosity.GE.verbosity_max) THEN
                  WRITE(6,'(a,3g14.6)') 'Translation displacement added:&
                      & u(1:3) = ', uTranslate(1:3)
                  WRITE(6,*)
          END IF
  END IF

  ! Clip atoms in unit cell applying periodic boundary conditions
  IF (clipAtom) THEN
          CALL Clip_Atoms(xp, im, at, 6, verbosity.GE.3)
          CALL Clip_Atoms(xp0, im0, at0, 6, verbosity.GE.3)
          clipAtom=.false.
  END IF

  ! Calculate homogeneous strain tensor
  IF (at_defined.AND.at0_defined) THEN
          CALL Mat3Inv(at0, inv_at0)
          epsi = 0.5d0*( MatMul( Transpose(inv_at0), MatMul( Transpose(at), &
                MatMul(at, inv_at0 ) ) ) - matId )
          IF (verbosity.GE.verbosity_max) THEN
                  WRITE(6,'(a)') 'Strain tensor corresponding to periodicity vectors'
                  DO i=1, 3
                      WRITE(6,'(a,3g14.6,a)') '      | ', epsi(i,1:3), ' |'
                  END DO
                  WRITE(6,*)
          END IF
  END IF


  ! Calculate standard deviation
  IF ( (verbosity.GE.verbosity_max).AND.out_displacement ) THEN
          uSolid(1:3)=0.d0 ; u2Solid(1:3)=0.d0
          DO n=1, im
                     uSolid(1:3) = uSolid(1:3) + uTotal(1:3,n)
                     u2Solid(1:3) = u2Solid(1:3) + uTotal(1:3,n)**2
          END DO
          IF (im.GT.0) uSolid(1:3) = uSolid(1:3)/dble(im)
          IF (im.GT.0) u2Solid(1:3) = u2Solid(1:3)/dble(im)
          WRITE(6,*)
          WRITE(6,'(a,3(g14.6,1x))') 'Average displacement (translation): U(1:3) = ', uSolid(:)
          WRITE(6,'(a,3(g14.6,1x))') 'Displacement standard deviation: dU(1:3) = ', &
                  sqrt( u2Solid(:) - uSolid(:)**2 )
          WRITE(6,'(a,g14.6)') '  dU = ', sqrt( Sum(u2Solid(:))/3.d0 - ( Sum(uSolid(:))/3.d0 )**2 )
          WRITE(6,'(a,g14.6)') 'Norm of the 3N vector (with translation)', sqrt( Sum( uTotal(1:3,1:im)**2 ) )
          dU2 = 0.d0
          DO n=1, im
                dU2 = dU2 + Sum( ( uTotal(1:3,n) - uSolid(1:3) )**2 )
          END DO
          WRITE(6,'(a,g14.6)') '                   (without translation)', sqrt( dU2 )
  END IF

  ! Calculate neighbours table for reference structure (if necessary)
  !IF (out_gradDisplacement.OR.out_pattern.OR.out_gradElasticDisplacement&
                !.OR.out_strain.OR.out_elasticStrain.OR.out_plasticStrain.OR.out_nye.OR.out_BurgersDensity) THEN
  IF (out_gradDisplacement.OR.out_strain.OR.out_plasticStrain) THEN
          CALL InitNeighbours(imm)
          ALLOCATE(nNeigh0(1:imm))               ; nNeigh0(:) = 0
          ALLOCATE(iNeigh0(1:max_nNeigh, 1:imm)) ; iNeigh0(:,:) = 0
          CALL BuildNeighbours(im0, xp0, at0, nNeigh0, iNeigh0, 6)
  END IF

  ! Calculate neighbours table for input structure (if necessary)
  IF (out_neighbours.OR.out_pattern.OR.out_gradElasticDisplacement&
                .OR.out_elasticStrain.OR.out_plasticStrain.OR.out_nye.OR.out_BurgersDensity) THEN
          CALL InitNeighbours(imm)
          ALLOCATE(nNeigh(1:imm))               ; nNeigh(:) = 0
          ALLOCATE(iNeigh(1:max_nNeigh, 1:imm)) ; iNeigh(:,:) = 0
          CALL BuildNeighbours(im, xp, at, nNeigh, iNeigh, 6)
  END IF

  ! Calculate gradient of total displacement
  IF (out_gradDisplacement.OR.out_strain.OR.out_plasticStrain) THEN
          CALL InitGradDisplacement(imm)
          CALL BuildGradDisplacement(im, xp0, uTotal, at0, at, nNeigh0, iNeigh0, 6)
  END IF

  ! Calculate atomic strain tensor (total tensor)
  IF (out_strain.OR.out_plasticStrain) THEN
          ALLOCATE(eAtomicStrain(1:6,1:imm)) ; eAtomicStrain(:,:)=0.d0
          CALL StrainFromDisplacement(im, gradDisplacement, eAtomicStrain, 6)
          IF (.NOT.(out_gradDisplacement)) CALL DestroyGradDisplacement()
  END IF

  ! Calculate gradient of elastic displacement
  IF (out_pattern.OR.out_gradElasticDisplacement.OR.out_elasticStrain.OR.out_plasticStrain.OR.out_nye.OR.out_BurgersDensity) THEN
          INQUIRE(file=patternFile, exist=ok)
          IF (.NOT.ok) THEN
                  WRITE(0,'(a)') 'You need to use variable "patternFile" to define a file'
                  WRITE(0,'(a)') '  where the crystallographic patterns defining the atom reference neighbourhoods are stored'
                  STOP '< Displacement >'
          END IF
          OPEN(file=patternFile, unit=55, status='old', action='read')
          CALL ReadPattern(55, 6)
          CLOSE(55)
          CALL InitGradElasticDisplacement(imm)
          CALL BuildGradElasticDisplacement(im, xp, at, nNeigh, iNeigh, 6)
  END IF

  ! Calculate atomic strain tensor (elastic part) 
  IF (out_elasticStrain.OR.out_plasticStrain) THEN
          ALLOCATE(eAtomicElasticStrain(1:6,1:imm)) ; eAtomicElasticStrain(:,:)=0.d0
          CALL StrainFromDisplacement(im, gradElasticDisplacement, eAtomicElasticStrain, 6)
          !IF (.NOT.(out_gradElasticDisplacement.OR.out_Nye.OR.out_BurgersDensity)) CALL DestroyGradElasticDisplacement()
  END IF

  ! Calculate atomic strain tensor (plastic part)
  IF (out_plasticStrain) THEN
          ALLOCATE(eAtomicPlasticStrain(1:6,1:imm)) ; eAtomicPlasticStrain(:,:)=0.d0
          eAtomicPlasticStrain(:,1:im) = eAtomicStrain(:,1:im) - eAtomicElasticStrain(:,1:im)
  END IF

  ! Calculate Nye tensor
  IF (out_nye.OR.out_BurgersDensity) THEN
          CALL InitNyeTensor(imm)
          !CALL BuildNyeTensor(im, xp0, gradElasticDisplacement, at0,  nNeigh, iNeigh, 6) ! TEST 2017/01/05
          CALL BuildNyeTensor(im, xp, gradElasticDisplacement, at,  nNeigh, iNeigh, 6)    ! TEST 2017/01/05
          !IF (.NOT.out_gradElasticDisplacement) CALL DestroyGradElasticDisplacement()
  END IF

  ! Properties printed on output
  nAux_real=0 ; nAux_int=0; n=0
  IF (out_neighbours)              nAux_int  = nAux_int  + 1
  IF (out_pattern)                 nAux_int  = nAux_int  + 1
  IF (out_displacement)            nAux_real = nAux_real + 3
  IF (out_gradDisplacement)        nAux_real = nAux_real + 9
  IF (out_gradElasticDisplacement) nAux_real = nAux_real + 9
  IF (out_Nye)                     nAux_real = nAux_real + 9
  IF (out_BurgersDensity)          nAux_real = nAux_real + 1
  IF (out_strain)                  nAux_real = nAux_real + 6
  IF (out_elasticStrain)           nAux_real = nAux_real + 8
  IF (out_plasticStrain)           nAux_real = nAux_real + 8

  IF ( (nAux_int.NE.0) .OR. (nAux_real.NE.0) ) THEN

          ALLOCATE(aux_title(1:nAux_real+nAux_int))
          IF (nAux_int.NE.0) THEN
                  ALLOCATE(aux_int(1:nAux_int,1:imm))
                  aux_int(:,:)=0
          END IF
          IF (nAux_real.NE.0) THEN
                  ALLOCATE(aux_real(1:nAux_real,1:imm)) 
                  aux_real(:,:) = 0.d0
          END IF

          IF (out_neighbours) THEN      ! Print number of neighbours for each atom
                  aux_int(n+1,:) = nNeigh(:)
                  aux_title(n+1)='Number of neighbours'
                  n=n+1
          END IF
          IF (out_pattern) THEN      ! Print pattern index for each atom
                  aux_int(n+1,1:im) = inpPattern(1:im)
                  aux_title(n+1)='Pattern index'
                  n=n+1
          END IF
          IF (out_displacement) THEN  ! Print atom displacement
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = uTotal(1:3,:)
                  aux_title(n+1)='Total displacement Ux'
                  aux_title(n+2)='                   Uy'
                  aux_title(n+3)='                   Uz (A)'
                  n=n+3
          END IF
          IF (out_gradDisplacement) THEN  ! Print gradient of total displacement
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = gradDisplacement(1,1:3,:)
                  aux_title(n+1)='Gradient of total displacement d Ux / dx'
                  aux_title(n+2)='                      d Ux / dy'
                  aux_title(n+3)='                      d Ux / dz'
                  n=n+3
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = gradDisplacement(2,1:3,:)
                  aux_title(n+1)='Gradient of total displacement d Uy / dx'
                  aux_title(n+2)='                      d Uy / dy'
                  aux_title(n+3)='                      d Uy / dz'
                  n=n+3
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = gradDisplacement(3,1:3,:)
                  aux_title(n+1)='Gradient of total displacement d Uz / dx'
                  aux_title(n+2)='                      d Uz / dy'
                  aux_title(n+3)='                      d Uz / dz'
                  n=n+3
          END IF
          IF (out_gradElasticDisplacement) THEN  ! Print gradient of eleastic displacement
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = gradElasticDisplacement(1,1:3,:)
                  aux_title(n+1)='Gradient of elastic displacement d Ux / dx'
                  aux_title(n+2)='                      d Ux / dy'
                  aux_title(n+3)='                      d Ux / dz'
                  n=n+3
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = gradElasticDisplacement(2,1:3,:)
                  aux_title(n+1)='Gradient of elastic displacement d Uy / dx'
                  aux_title(n+2)='                      d Uy / dy'
                  aux_title(n+3)='                      d Uy / dz'
                  n=n+3
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = gradElasticDisplacement(3,1:3,:)
                  aux_title(n+1)='Gradient of elastic displacement d Uz / dx'
                  aux_title(n+2)='                      d Uz / dy'
                  aux_title(n+3)='                      d Uz / dz'
                  n=n+3
          END IF
          IF (out_Nye) THEN  ! Print Nye tensor
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = nyeAtomic(1,1:3,:)
                  aux_title(n+1)='Nye tensor a_x,x '
                  aux_title(n+2)='           a_x,y'
                  aux_title(n+3)='           a_x,z'
                  n=n+3
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = nyeAtomic(2,1:3,:)
                  aux_title(n+1)='Nye tensor a_y,x '
                  aux_title(n+2)='           a_y,y'
                  aux_title(n+3)='           a_y,z'
                  n=n+3
                  aux_real(n-nAux_int+1:n-nAux_int+3,:) = nyeAtomic(3,1:3,:)
                  aux_title(n+1)='Nye tensor a_z,x '
                  aux_title(n+2)='           a_z,y'
                  aux_title(n+3)='           a_z,z'
                  n=n+3
          END IF
          IF (out_burgersDensity) THEN  ! Print Burgers density defined by bNye(:) and lNye(:)
                  DO i=1, im
                     aux_real(n-nAux_int+1:n-nAux_int+1,i) = Sum( bNye(1:3) * MatMul( nyeAtomic(1:3,1:3,i), lNye(1:3) ) )
                  END DO
                  aux_title(n+1)='Burgers density '
                  n=n+1
          END IF
          IF (out_strain) THEN        ! Print atomic strain (6 components)
                  aux_real(n-nAux_int+1:n-nAux_int+6,:)=eAtomicStrain(1:6,:)
                  aux_title(n+1) = 'Total strain Exx'
                  aux_title(n+2) = '       Eyy'
                  aux_title(n+3) = '       Ezz'
                  aux_title(n+4) = '     2*Eyz'
                  aux_title(n+5) = '     2*Exz'
                  aux_title(n+6) = '     2*Exy'
                  n=n+6
          END IF
          IF (out_elasticStrain) THEN        ! Print atomic strain (6 components)
                  aux_real(n-nAux_int+1:n-nAux_int+6,:)=eAtomicElasticStrain(1:6,:)
                  aux_real(n-nAux_int+7,:) = eAtomicElasticStrain(1,:) + eAtomicElasticStrain(2,:) + eAtomicElasticStrain(3,:)
                  aux_real(n-nAux_int+8,:) = Sqrt( 0.5d0*( ( eAtomicElasticStrain(1,:) - eAtomicElasticStrain(2,:) )**2 &
                         + ( eAtomicElasticStrain(2,:) - eAtomicElasticStrain(3,:) )**2 & 
                         + ( eAtomicElasticStrain(3,:) - eAtomicElasticStrain(2,:) )**2 & 
                         + 6.d0*( eAtomicElasticStrain(4,:)**2 + eAtomicElasticStrain(5,:)**2 + eAtomicElasticStrain(6,:)**2 ) ) )
                  aux_title(n+1) = 'Elastic strain Exx'
                  aux_title(n+2) = '       Eyy'
                  aux_title(n+3) = '       Ezz'
                  aux_title(n+4) = '     2*Eyz'
                  aux_title(n+5) = '     2*Exz'
                  aux_title(n+6) = '     2*Exy'
                  aux_title(n+7) = '     Tr(Eij)'
                  aux_title(n+8) = '     |dev(Eij)| * sqrt(3/2).'
                  n=n+8
          END IF
          IF (out_plasticStrain) THEN        ! Print atomic strain (6 components)
                  aux_real(n-nAux_int+1:n-nAux_int+6,:)=eAtomicPlasticStrain(1:6,:)
                  aux_real(n-nAux_int+1:n-nAux_int+6,:)=eAtomicPlasticStrain(1:6,:)
                  aux_real(n-nAux_int+7,:) = eAtomicPlasticStrain(1,:) + eAtomicPlasticStrain(2,:) + eAtomicPlasticStrain(3,:)
                  aux_real(n-nAux_int+8,:) = Sqrt( 0.5d0*( ( eAtomicPlasticStrain(1,:) - eAtomicPlasticStrain(2,:) )**2 &
                         + ( eAtomicPlasticStrain(2,:) - eAtomicPlasticStrain(3,:) )**2 & 
                         + ( eAtomicPlasticStrain(3,:) - eAtomicPlasticStrain(2,:) )**2 & 
                         + 6.d0*( eAtomicPlasticStrain(4,:)**2 + eAtomicPlasticStrain(5,:)**2 + eAtomicPlasticStrain(6,:)**2 ) ) )
                  aux_title(n+1) = 'Plastic strain Exx'
                  aux_title(n+2) = '       Eyy'
                  aux_title(n+3) = '       Ezz'
                  aux_title(n+4) = '     2*Eyz'
                  aux_title(n+5) = '     2*Exz'
                  aux_title(n+6) = '     2*Exy'
                  aux_title(n+7) = '     Tr(Eij)'
                  aux_title(n+8) = '     |dev(Eij)| * sqrt(3/2).'
                  n=n+8
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
          IF (outXyz) WRITE(6,'(2a)') 'Write output structure with Xyz format in file ', outFile
          IF (outOnlyAtoms) WRITE(6,'(2a)') 'Write one line per atom in file ', outFile
          IF (outCfg) WRITE(6,'(2a)') 'Write output structure with Cfg format in file ', outFile
          IF (outLammps) WRITE(6,'(2a)') 'Write output structure with Lammps dump format in file ', outFile
          IF (outPoscar) WRITE(6,'(2a)') 'Write output structure with Poscar dump format in file ', outFile
          IF (initial) THEN
                  WRITE(6,'(a)') '  Initial atom coordinates are kept (reference structure)'
                  IF (at_defined) &
                        WRITE(6,'(a)') '  Initial periodicity vectors are kept'
          ELSE
                  WRITE(6,'(a)') '  Final atom coordinates are kept (input structure)'
                  IF (at_defined) &
                        WRITE(6,'(a)') '  Final periodicity vectors are kept'
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

  IF ( outXyz.OR.outOnlyAtoms.OR.outCfg.OR.outGin.OR.outLisa.OR.outSiesta.OR.outNDM.OR.outLammps.OR.outPoscar ) THEN

          IF (outXyz) fileType='xyz'
          IF (outOnlyAtoms) fileType='onlyAtoms'
          IF (outCfg) fileType='cfg'
          IF (outGin) WRITE(0,'(a)') "Gin type not supported for output structure file"
          IF (outLisa) WRITE(0,'(a)') "Lisa type not supported for output structure file"
          IF (outSiesta) WRITE(0,'(a)') "Siesta type not supported for output structure file"
          IF (outNDM) WRITE(0,'(a)') "Binary NDM type not supported for output structure file"
          IF (outLammps) fileType='lammps'
          IF (outPoscar) fileType='poscar'
          IF (initial) THEN
                  IF ( (nAux_int.GT.0).AND.(nAux_real.GT.0) ) THEN
                          CALL WriteStructure(xp0, iTyp0, im0, at0, outFile, fileType, &
                                nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                aux_title=aux_title(:))
                  ELSE IF (nAux_int.GT.0) THEN
                          CALL WriteStructure(xp0, iTyp0, im0, at0, outFile, fileType, &
                                nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                aux_title=aux_title(:))
                  ELSE IF (nAux_real.GT.0) THEN
                          CALL WriteStructure(xp0, iTyp0, im0, at0, outFile, fileType, &
                                nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                aux_title=aux_title(:))
                  ELSE
                          CALL WriteStructure(xp0, iTyp0, im0, at0, outFile, fileType)
                  END IF
          ELSE
                  IF ( (nAux_int.GT.0).AND.(nAux_real.GT.0) ) THEN
                          CALL WriteStructure(xp, iTyp, im, at, outFile, fileType, &
                                nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                aux_title=aux_title(:))
                  ELSE IF (nAux_int.GT.0) THEN
                          CALL WriteStructure(xp, iTyp, im, at, outFile, fileType, &
                                nAux_int=nAux_int, aux_int=aux_int(:,:), &
                                aux_title=aux_title(:))
                  ELSE IF (nAux_real.GT.0) THEN
                          CALL WriteStructure(xp, iTyp, im, at, outFile, fileType, &
                                nAux_real=nAux_real, aux_real=aux_real(:,:), &
                                aux_title=aux_title(:))
                  ELSE
                          CALL WriteStructure(xp, iTyp, im, at, outFile, fileType)
                  END IF
          END IF

  END IF


END PROGRAM Displacement
