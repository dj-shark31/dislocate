MODULE EnergyDipoleModule

  INTEGER, parameter, private :: verbosity_max=2

CONTAINS

  SUBROUTINE ElasticEnergyDipole(Etotal, Serr_DDipole, out)
    ! Calculate elastic energy for a distribution of dislocation dipoles

    USE Babel_Data
    USE DDipoleModule
    USE PeriodModule
    USE Math

    IMPLICIT NONE
    REAL(kind(0.d0)), intent(out) :: Etotal        ! Elastic energy
    ! Correction to stress tensor due to conditional convergence for
    ! dislocations and line-force couples
    REAL(kind(0.d0)), dimension(1:6), intent(in) :: Serr_DDipole
    INTEGER, intent(in) :: out

    INTEGER :: n, m, ix, iy, iz
    REAL(kind(0.d0)) :: E0, Edipole, Eerr, K0
    REAL(kind(0.d0)) :: Einter_dipole_dipole_pr, Einter_dipole_dipole_im, &
        Eel_strain, Einter_dipole_strain
    REAL(kind(0.d0)) :: volume, length, dR2_11, dR2_12, dR2_21, dR2_22
    REAL(kind(0.d0)), dimension(1:3) :: R, dA
    REAL(kind(0.d0)), dimension(1:3,1:3) :: Sigma, sHomogeneous
    REAL(kind(0.d0)), dimension(1:6) :: sVoigt

    ! Initialization
    Etotal=0.d0

    ! Check for dislocation dipoles
    IF (l_dDipole) THEN
            IF (.NOT.one_DDipole_system) THEN
                    IF (verbosity.GE.verbosity_max) THEN
                            WRITE(out,'(a)') 'Elastic energy cannot be calculated as all&
                                & dislocation dipoles do not have the same line direction'
                            WRITE(out,*)
                            WRITE(out,'(a)') '==========================='
                            WRITE(out,*)
                    END IF
                    RETURN
            END IF
            DO n=1, nDDipole
               DO m=n+1, nDDipole
                  dR2_11 = Sum( ( c1DDipole(:,n) - c1DDipole(:,m) )**2 )
                  dR2_12 = Sum( ( c1DDipole(:,n) - c2DDipole(:,m) )**2 )
                  dR2_21 = Sum( ( c2DDipole(:,n) - c1DDipole(:,m) )**2 )
                  dR2_22 = Sum( ( c2DDipole(:,n) - c2DDipole(:,m) )**2 )
                  IF ( ( dR2_11.LE.distance_zero2 ).OR.( dR2_12.LE.distance_zero2 ) &
                        .OR.( dR2_21.LE.distance_zero2 ).OR.( dR2_22.LE.distance_zero2 ) ) THEN
                          IF (verbosity.GE.verbosity_max) THEN
                                  WRITE(out,'(a,2(i0,a))') 'Elastic energy cannot be calculated as some&
                                      & dislocations of dipoles ', n, ' and ', m,' share the same position'
                                  WRITE(out,*)
                                  WRITE(out,'(a)') '==========================='
                                  WRITE(out,*)
                          END IF
                          RETURN
                  END IF
               END DO
            END DO
    END IF

    ! Check for homogeneous strain
    IF (strain.AND..NOT.at_defined) THEN
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a)') 'Elastic energy cannot be calculated as&
                        & homogenenous strain is applied and no unit cell is&
                        & defined (infinite volume)'
                    WRITE(out,*)
                    WRITE(out,'(a)') '==========================='
                    WRITE(out,*)
            END IF
            RETURN
    END IF

    ! Unit length for line defects
    IF (.NOT.at_defined) THEN
            length = 1.d0
    ELSE IF (l_DDipole) THEN
            IF (Colinear(at(:,3),lDDipole(:,1))) THEN
                    length = Sqrt( Sum( at(1:3,3)**2 ) )
            ELSE IF (Colinear(at(:,2),lDDipole(:,1))) THEN
                    length = Sqrt( Sum( at(1:3,2)**2 ) )
            ELSE IF (Colinear(at(:,1),lDDipole(:,1))) THEN
                    length = Sqrt( Sum( at(1:3,1)**2 ) )
            ELSE
                    length = 1.d0
            END IF
    ELSE
            length = 1.d0
    END IF
            
    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,'(a,g14.6)') ' Elastic energy calculation with core cutoff radius  rc = ', rc
            WRITE(out,'(a,g14.6)') '   factor used to convert [Stress]*[Distance]^3 in [Energy]: &
                & factorE = ', factorE
            WRITE(out,'(a,g14.6,a)') '   ( 1 GPa * 1 A^3 = ', 1.d9*1.d-30/1.602176530000000045d-19,' eV )'
            WRITE(out,'(a,g14.6)') '  unit length for line defects = ', length
            WRITE(out,*)
    END IF

    ! Initialisation
    Etotal = 0.d0
    Einter_dipole_dipole_pr = 0.d0
    Einter_dipole_dipole_im = 0.d0
    Eel_strain = 0.d0
    Einter_dipole_strain = 0.d0

    ! Elastic energy due to homogeneous strain
    IF (strain) THEN
            volume=abs(MatDet(at))
            ! Homogeneous strain
            sVoigt(1:6) = MatMul( CVoigt, (/ eStrain(1,1), eStrain(2,2), eStrain(3,3), &
                eStrain(2,3)+eStrain(3,2), eStrain(1,3)+eStrain(3,1), eStrain(1,2)+eStrain(2,1) /) )
            sHomogeneous = Reshape( (/ sVoigt(1), sVoigt(6), sVoigt(5), &
                                        sVoigt(6), sVoigt(2), sVoigt(4), &
                                        sVoigt(5), sVoigt(4), sVoigt(3) /), (/3,3/) )
            Eel_strain = factorE*0.5d0*volume * Sum( sHomogeneous(1:3,1:3)*eStrain(1:3,1:3) )
    END IF

    IF (l_DDipole) THEN
            ! Elastic energy of dislocation dipoles located in primitive unit cell
            dislo_loop1: DO n=1, nDDipole
               ! Indexes of dislocations forming a dipole
               Edipole = length*factorE*DDipole_Self_Energy(n, rc)
               Etotal = Etotal + Edipole
               IF (verbosity.GE.verbosity_max) THEN
                       WRITE(out,'(a,i0)') '  Dipole ', n
                       WRITE(out,'(2(a,g16.8),a)') '    Edipole = E0 - K0 ln(rc) = &
                                &2 A(theta) + K0 ln(|d|/rc) = ', Edipole, ' = ', &
                                Edipole/length, ' per unit length'
                       WRITE(out,*)
               END IF
            END DO dislo_loop1
        
            ! Elastic interaction energy between dislocations belonging to different dipoles
            dislo_loop2: DO n=1, nDDipole
                  DO m=n+1, nDDipole
                     Edipole = length*factorE*DDipole_Interaction_Energy(n,m)
                     IF (verbosity.GE.verbosity_debug) THEN
                             WRITE(out,'(2(a,i0),2(a,g16.8),a)') '  Interaction energy between&
                                      & dislocation dipoles ', n, ' and ', &
                                     m, ':  Einter = ', Edipole, ' = ', Edipole/length, &
                                     ' per unit length'
                     END IF
                     Einter_dipole_dipole_pr = Einter_dipole_dipole_pr + Edipole
                  END DO
            END DO dislo_loop2
    
            ! Interaction with image dislocations
            images_loop: DO ix=-nxImages, nxImages      ! Loop on all image cells
              DO iy=-nyImages, nyImages
                 DO iz=-nzImages, nzImages
                   IF ( (ix.EQ.0).AND.(iy.EQ.0).AND.(iz.EQ.0) ) Cycle
                   R(1:3) = dble(ix)*at(1:3,1) + dble(iy)*at(1:3,2) &
                        + dble(iz)*at(1:3,3)
                   DO n=1, nDDipole
                      DO m=1, nDDipole
                         Einter_dipole_dipole_im = Einter_dipole_dipole_im &
                                + DDipole_Interaction_Energy(n, m, shift=R)
                      END DO
                   END DO
            END DO ; END DO ; END DO images_loop
            Einter_dipole_dipole_im = 0.5d0*length*factorE*Einter_dipole_dipole_im

            ! Correction to interaction with image dislocations
            !  and interaction with homogeneous strain
            sigma(1,1) = Serr_DDipole(1) ; sigma(2,2) = Serr_DDipole(2) ; sigma(3,3) = Serr_DDipole(3)
            sigma(2,3) = Serr_DDipole(4) ; sigma(1,3) = Serr_DDipole(5) ; sigma(1,2) = Serr_DDipole(6)
            sigma(3,2) = Serr_DDipole(4) ; sigma(3,1) = Serr_DDipole(5) ; sigma(2,1) = Serr_DDipole(6)
            Eerr = 0.d0
            dislo_loop3: DO n=1, nDDipole
               dA(:) = -CrossProduct( lDDipole(:,n), c1DDipole(:,n)-c2DDipole(:,n) )
               Eerr = Eerr + Sum( bDDipole(:,n)*MatMul( Sigma(:,:), dA ) )
               IF (strain) Einter_dipole_strain = Einter_dipole_strain &
                         + Sum( bDDipole(:,n)*MatMul( sHomogeneous(:,:), dA ) )
            END DO dislo_loop3
            Einter_dipole_dipole_im = Einter_dipole_dipole_im + 0.5d0*length*factorE*Eerr
            Einter_dipole_strain = length*factorE*Einter_dipole_strain
    END IF ! Dislocation

    ! Add different contributions
    IF (.NOT.l_DDipole) THEN
            K0=0.d0
    ELSE
            K0 = length*factorE*2.d0*Sum( DDipole_prelog_energy_factor(1:nDDipole) )
    END IF
    Etotal = Etotal + Einter_dipole_dipole_pr + Einter_dipole_dipole_im + Einter_dipole_strain + Eel_strain
    E0 = Etotal + K0*log(rc)
    
    IF (verbosity.GE.verbosity_max) THEN
            IF (strain) THEN
                    WRITE(out,'(2(a,g16.8),a)') '  Homogeneous strain contribution: &
                        &                    E0 = ', Eel_strain, " = ", Eel_strain/length, " per unit length"
                    IF (l_DDipole) WRITE(out,'(2(a,g16.8),a)') &
                        '                 interaction with dislocation dipoles: E0 = ', &
                        Einter_dipole_strain, ' = ', Einter_dipole_strain/length, ' per unit length'
                    WRITE(out,*)
            END IF
            IF (l_DDipole.AND.(nDDipole.GE.2)) THEN
                    WRITE(out,'(a)') '  Interaction between line-defects contained in primitive unit cell'
                    WRITE(out,'(2(a,g16.8),a)') &
                       '         dislocation dipoles - dislocation dipoles:    E0 = ', &
                       Einter_dipole_dipole_pr, ' = ', Einter_dipole_dipole_pr/length, ' per unit length'
                    WRITE(out,*)
            END IF
            IF (l_DDipole.AND.(xImages.OR.yImages.OR.zImages)) THEN
                    WRITE(out,'(a)') '  Dislocation interaction with image line-defects'
                    WRITE(out,'(2(a,g16.8),a)') '                                image dislocations:    E0 = ', &
                       Einter_dipole_dipole_im, ' = ', Einter_dipole_dipole_im/length, ' per unit length'
                    WRITE(out,*)
            END IF
            WRITE(out,'(a)') '--------------'
            WRITE(out,*)
            WRITE(out,'(2(a,g16.8),a)') ' Total elastic energy :    E = E0 - K0 ln(rc) - K2 / rc^2 = ', Etotal, ' = ', Etotal/length, ' per unit length'
            WRITE(out,'(2(a,g16.8),a)') '       with E0 = ', E0, ' = ', E0/length, ' per unit length'
            WRITE(out,'(2(a,g16.8),a)') '            K0 = ', K0, ' = ', K0/length, ' per unit length'
            WRITE(out,'(2(a,g16.8),a)') '            K2 = ', 0., ' = ', 0.,        ' per unit length'
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
    END IF

  END SUBROUTINE ElasticEnergyDipole

END MODULE EnergyDipoleModule
