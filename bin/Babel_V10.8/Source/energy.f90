MODULE EnergyModule

  INTEGER, parameter, private :: verbosity_max=2

CONTAINS

  SUBROUTINE ElasticEnergy(Etotal, Serr_dislo, Serr_lineCouple, out)
    ! Calculate elastic energy
    ! This only works for dislocations if they all have
    ! the same line direction (only 1 system), the total Burgers
    ! vector is null and dislocations can be grouped in dipoles.

    USE Babel_Data
    USE Rearrange
    USE Disloc_elasticity_ani
    USE LineCouple_elasticity_ani
    USE LineForce_elasticity_ani
    USE PeriodModule
    USE Math

    IMPLICIT NONE
    REAL(kind(0.d0)), intent(out) :: Etotal        ! Elastic energy
    ! Correction to stress tensor due to conditional convergence for
    ! dislocations and line-force couples
    REAL(kind(0.d0)), dimension(1:6), intent(in) :: Serr_dislo, Serr_lineCouple
    INTEGER, intent(in) :: out

    INTEGER :: n, m, n1, n2, n3, n4, ix, iy, iz, nd_temp, nlc_temp
    REAL(kind(0.d0)) :: E0, E0dipole, Eerr, K0, K2
    REAL(kind(0.d0)) :: Einter_dipole_dipole_pr, Einter_dipole_couple_pr, Einter_couple_couple_pr,&
        Einter_dipole_dipole_im, Einter_couple_dipole_im, Einter_couple_couple_im, &
        Eel_strain, Einter_dipole_strain, Einter_couple_strain
    REAL(kind(0.d0)) :: volume, length, d, A, dR2
    REAL(kind(0.d0)), dimension(1:3) :: R, R0, dA,dR
    REAL(kind(0.d0)), dimension(1:3,1:3) :: Sigma, sHomogeneous, epsi
    REAL(kind(0.d0)), dimension(1:6) :: sVoigt

    IF (lineForce) THEN
            IF (verbosity.GE.verbosity_max) THEN
                    WRITE(out,'(a)') 'Elastic energy not implemented with&
                        & line-forces' 
                    WRITE(out,*)
                    WRITE(out,'(a)') '==========================='
                    WRITE(out,*)
            END IF
            RETURN
    END IF

    ! Initialization
    Etotal=0.d0

    ! Check for dislocations
    IF (l_dislo) THEN
            IF (nds.GT.1) THEN
                    IF (verbosity.GE.verbosity_max) THEN
                            WRITE(out,'(a)') 'Elastic energy cannot be calculated as all&
                                & dislocations do not have the same line direction'
                            WRITE(out,*)
                            WRITE(out,'(a)') '==========================='
                            WRITE(out,*)
                    END IF
                    RETURN
            ELSEIF (.NOT.bClosedSystem(1)) THEN
                    IF (verbosity.GE.verbosity_max) THEN
                            WRITE(out,'(a)') 'Elastic energy cannot be calculated as&
                                & the total Burgers vector is not null'
                            WRITE(out,*)
                            WRITE(out,'(a)') '==========================='
                            WRITE(out,*)
                    END IF
                    RETURN
            ELSEIF(.NOT.bDipoleSystem(1)) THEN
                    IF (verbosity.GE.verbosity_max) THEN
                            WRITE(out,'(a)') 'Elastic energy cannot be calculated as&
                                & dislocations could not be grouped in dipoles'
                            WRITE(out,*)
                            WRITE(out,'(a)') '==========================='
                            WRITE(out,*)
                    END IF
                    RETURN
            END IF
            DO n=1, nd
               DO m=n+1, nd
                  dR2 = Sum( ( cDislo(:,n) - cDislo(:,m) )**2 )
                  IF ( dR2.LE.distance_zero2 ) THEN
                          IF (verbosity.GE.verbosity_max) THEN
                                  WRITE(out,'(a,2(i0,a))') 'Elastic energy cannot be calculated as&
                                      & dislocations ', n, ' and ', m,' share the same position'
                                  WRITE(out,*)
                                  WRITE(out,'(a)') '==========================='
                                  WRITE(out,*)
                          END IF
                          RETURN
                  END IF
               END DO
            END DO
    END IF

    ! Check for line-force couples
    IF (l_lineCouple) THEN
            IF (nlcs.GT.1) THEN
                    IF (verbosity.GE.verbosity_max) THEN
                            WRITE(out,'(a)') 'Elastic energy cannot be calculated as all&
                                & line-force couples do not have the same line direction'
                            WRITE(out,*)
                            WRITE(out,'(a)') '==========================='
                            WRITE(out,*)
                    END IF
                    RETURN
            END IF
    END IF

    ! Check for dislocations and line-force couples combination
    IF (l_dislo.AND.l_lineCouple) THEN
            IF (.NOT.Colinear(lDisloSystem(:,1),lLineCoupleSystem(:,1))) THEN
                    IF (verbosity.GE.verbosity_max) THEN
                            WRITE(out,'(a)') 'Elastic energy cannot be calculated as dislocations&
                                & and line-force couples do not share the same line direction'
                            WRITE(out,*)
                            WRITE(out,'(a)') '==========================='
                            WRITE(out,*)
                    END IF
                    RETURN
            END IF
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
    ELSE IF (l_dislo) THEN
            IF (Colinear(at(:,3),lDisloSystem(:,1))) THEN
                    length = Sqrt( Sum( at(1:3,3)**2 ) )
            ELSE IF (Colinear(at(:,2),lDisloSystem(:,1))) THEN
                    length = Sqrt( Sum( at(1:3,2)**2 ) )
            ELSE IF (Colinear(at(:,1),lDisloSystem(:,1))) THEN
                    length = Sqrt( Sum( at(1:3,1)**2 ) )
            ELSE
                    length = 1.d0
            END IF
    ELSE IF (l_lineCouple) THEN
            IF (Colinear(at(:,3),lLineCoupleSystem(:,1))) THEN
                    length = Sqrt( Sum( at(1:3,3)**2 ) )
            ELSE IF (Colinear(at(:,2),lLineCoupleSystem(:,1))) THEN
                    length = Sqrt( Sum( at(1:3,2)**2 ) )
            ELSE IF (Colinear(at(:,1),lLineCoupleSystem(:,1))) THEN
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
    E0 = 0.d0
    Einter_dipole_dipole_pr=0.d0
    Einter_dipole_couple_pr=0.d0
    Einter_couple_couple_pr=0.d0
    Einter_dipole_dipole_im=0.d0
    Einter_couple_dipole_im=0.d0
    Einter_couple_couple_im=0.d0
    Eel_strain = 0.d0
    Einter_couple_strain=0.d0
    Einter_dipole_strain=0.d0

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

    IF (l_dislo) THEN
            ! Elastic energy of dislocation dipoles located in primitive unit cell
            dislo_loop1: DO n=1, nDisloSystem(1), 2
               ! Indexes of dislocations forming a dipole
               n1=iDisloSystem(n,1)
               n2=iDisloSystem(n+1,1)
               E0dipole = length*factorE*( Dislo_Dislo_Interaction_ani(n1,n2) &
                       - dislo_prelog_energy_factor(n1)*log(rc) &
                       - dislo_prelog_energy_factor(n2)*log(rc) &
                       + dislo_angular_energy(n1) + dislo_angular_energy(n2) )
               Etotal = Etotal + E0dipole
               IF (verbosity.GE.verbosity_max) THEN
                       K0 = length*factorE*( dislo_prelog_energy_factor(n1) + dislo_prelog_energy_factor(n2) )
                       dR(:) = cDislo(:,n1) - cDislo(:,n2) &
                               - Sum( ( cDislo(:,n1) - cDislo(:,n2) )*lDislo(:,n1) )*lDislo(:,n1)
                       d = Sqrt( Sum( dR(:)**2 ) )
                       A = 0.5d0*( E0dipole - K0*log(d/rc) )
                       WRITE(out,'(2(a,i0))') '  Dipole formed by dislocations ', n1, ' and ', n2
                       WRITE(out,'(2(a,g16.8),a)') '    Edipole = E0 - K0 ln(rc) = 2 A(theta) + K0 ln(|d|/rc) = ', & 
                               E0dipole, ' = ', E0dipole/length, ' per unit length'
                       WRITE(out,'(2(a,g16.8),a)') '       with E0 = ', E0dipole + K0*log(rc), ' = ', (E0dipole + K0*log(rc))/length, ' per unit length'
                       WRITE(out,'(2(a,g16.8),a)') '       with K0 = ', K0,       ' = ', K0/length,       ' per unit length'
                       WRITE(out,'(2(a,g16.8),a)') '       with A  = ', A,        ' = ', A/length,        ' per unit length'
                       WRITE(out,'(a,g16.8)')      '       with d  = ', d
                       WRITE(out,*)
               END IF
            END DO dislo_loop1
        
            ! Elastic interaction energy between dislocations belonging to different dipoles
            dislo_loop2: DO n=1, nDisloSystem(1), 2
                  n1 = iDisloSystem(n,1)
                  n2 = iDisloSystem(n+1,1)
                  DO m=n+2, nDisloSystem(1), 2
                     n3 = iDisloSystem(m,1)
                     n4 = iDisloSystem(m+1,1)
                     E0dipole = length*factorE*Dipole_Dipole_Interaction_ani(n1,n2,n3,n4)
                     IF (verbosity.GE.verbosity_debug) THEN
                             WRITE(out,'(2(a,i0),a,g16.8)') '  Interaction energy between&
                                & dislo ', n1, ' and ', n3, ': ', &
                                length*factorE*Dislo_Dislo_Interaction_ani(n1, n3)
                             WRITE(out,'(2(a,i0),a,g16.8)') '  Interaction energy between&
                                & dislo ', n1, ' and ', n4, ': ', &
                                length*factorE*Dislo_Dislo_Interaction_ani(n1, n4)
                             WRITE(out,'(2(a,i0),a,g16.8)') '  Interaction energy between&
                                & dislo ', n2, ' and ', n3, ': ', &
                                length*factorE*Dislo_Dislo_Interaction_ani(n2, n3)
                             WRITE(out,'(2(a,i0),a,g16.8)') '  Interaction energy between&
                                & dislo ', n2, ' and ', n4, ': ', &
                                length*factorE*Dislo_Dislo_Interaction_ani(n2, n4)
                             WRITE(out,'(4(a,i0),2(a,g16.8),a)') '  Interaction energy between&
                                      & dislocation dipoles ', n1,'-',n2, ' and ', &
                                     n3,'-',n4, ':  Einter = ', E0dipole, ' = ', E0dipole/length, &
                                     ' per unit length'
                     END IF
                     Einter_dipole_dipole_pr = Einter_dipole_dipole_pr + E0dipole
                  END DO
            END DO dislo_loop2
    
            ! Interaction with image dislocations
            Einter_dipole_dipole_im = 0.d0
            images_loop: DO ix=-nxImages, nxImages      ! Loop on all image cells
              DO iy=-nyImages, nyImages
                 DO iz=-nzImages, nzImages
                   IF ( (ix.EQ.0).AND.(iy.EQ.0).AND.(iz.EQ.0) ) Cycle
                   R(1:3) = dble(ix)*at(1:3,1) + dble(iy)*at(1:3,2) &
                        + dble(iz)*at(1:3,3)
                   DO n=1, nDisloSystem(1), 2
                      n1 = iDisloSystem(n,1)
                      n2 = iDisloSystem(n+1,1)
                      DO m=1, nDisloSystem(1), 2
                         n3 = iDisloSystem(m,1)
                         n4 = iDisloSystem(m+1,1)
                         Einter_dipole_dipole_im = Einter_dipole_dipole_im &
                                + Dipole_Dipole_Interaction_ani(n1,n2, n3,n4, shift=R)
                      END DO
                   END DO
            END DO ; END DO ; END DO images_loop
            Einter_dipole_dipole_im = 0.5d0*length*factorE*Einter_dipole_dipole_im

            ! Correction to interaction with image dislocations
            !  and interaction with homogeneous strain
            sigma(1,1) = Serr_dislo(1) ; sigma(2,2) = Serr_dislo(2) ; sigma(3,3) = Serr_dislo(3)
            sigma(2,3) = Serr_dislo(4) ; sigma(1,3) = Serr_dislo(5) ; sigma(1,2) = Serr_dislo(6)
            sigma(3,2) = Serr_dislo(4) ; sigma(3,1) = Serr_dislo(5) ; sigma(2,1) = Serr_dislo(6)
            Eerr = 0.d0
            dislo_loop3: DO n=1, nDisloSystem(1), 2
               ! Indexes of dislocations forming a dipole
               n1=iDisloSystem(n,1)
               n2=iDisloSystem(n+1,1)
               dA(:) = -CrossProduct( lDislo(:,n1), cDislo(:,n1)-cDislo(:,n2) )
               Eerr = Eerr + Sum( bDislo(:,n1)*MatMul( Sigma(:,:), dA ) )
               IF (strain) Einter_dipole_strain = Einter_dipole_strain &
                         + Sum( bDislo(:,n1)*MatMul( sHomogeneous(:,:), dA ) )
            END DO dislo_loop3
            Einter_dipole_dipole_im = Einter_dipole_dipole_im + 0.5d0*length*factorE*Eerr
            Einter_dipole_strain = length*factorE*Einter_dipole_strain
    END IF ! Dislocation

    ! Line-couples 
    couple_loop1: DO n1=1, nlc
       IF (verbosity.GE.verbosity_max) THEN
               K2 = length*factorE*lineCouple_energy_factor(n1)
               WRITE(out,'(a,i0)') '  Line-force couple ', n1
               WRITE(out,'(2(a,g16.8),a)') '    E = -K2 / rc^2 = ', -K2/rc**2, ' = ', -K2/(rc**2*length), ' per unit length'
               WRITE(out,'(2(a,g16.8),a)') '       with K2  = ', K2, ' = ',  K2/length, ' per unit length'
               WRITE(out,*)
       END IF
       R0(1:3) = cLineCouple(1:3,n1)

       !-------------------------------------------------------------------
       ! Interaction with other line defects: calculate the stress created by
       ! these line defects 
       !   dislocations contained in primitive unit cell (*1)
       sVoigt(:) = 0.d0
       DO n2=1, nd
        ! Check if dislocation and line-force couples have the same location
         IF ( dislo_lineCouple(n2,n1) ) Cycle
         sVoigt(:) = sVoigt(:) + dislo_stress_ani(R0,n2)
       END DO
       Einter_dipole_couple_pr = Einter_dipole_couple_pr + factorE*length*Stress_LineCouple_Interaction(sVoigt,n1)

       !   line-force couples contained in primitive unit cell (*1/2)
       sVoigt(:) = 0.d0
       DO n2=n1+1, nlc
          sVoigt(:) = sVoigt(:) + lineCouple_stress_ani(R0,n2)
       END DO
       Einter_couple_couple_pr = Einter_couple_couple_pr + factorE*length*Stress_LineCouple_Interaction(sVoigt,n1)

       !   image dislocations (*1)
       nd_temp=nd ; nlc_temp=nlc
       nd=nd_temp ; nlc=0
       sVoigt(:) = Array_stress(R0,only_images=.true.) + Serr_dislo(:)
       Einter_couple_dipole_im = Einter_couple_dipole_im + factorE*length*Stress_LineCouple_Interaction(sVoigt,n1)

       !   image line-force couples (*1/2)
       nd=0 ; nlc=nlc_temp
       sVoigt(:) = 0.5d0*( Array_stress(R0,only_images=.true.) + Serr_lineCouple(:) )
       nd=nd_temp ; nlc=nlc_temp
       Einter_couple_couple_im = Einter_couple_couple_im + factorE*length*Stress_LineCouple_Interaction(sVoigt,n1)

       !   homogeneous strain
       epsi(:,:) = MatMul( rotLineCouple(:,:,n1), MatMul( eStrain(:,:), inv_rotLineCouple(:,:,n1) ) )
       Einter_couple_strain = Einter_couple_strain &
                 -factorE*length*( epsi(1,1)*mxLineCouple(n1)+epsi(2,2)*myLineCouple(n1) )
    END DO couple_loop1


    ! Add different contributions
    IF ((.NOT.l_dislo).OR.(nd.LE.0)) THEN
            K0=0.d0
    ELSE
            K0 = length*factorE*Sum( dislo_prelog_energy_factor(1:nd) )
    END IF
    IF ((.NOT.l_lineCouple).AND.(nlc.LE.0)) THEN
            K2=0.d0
    ELSE
            K2 = length*factorE*Sum( lineCouple_energy_factor(1:nlc) )
    END IF
    Etotal = Etotal + Einter_dipole_dipole_pr + Einter_dipole_dipole_im + Einter_dipole_strain &
            + Einter_dipole_couple_pr + Einter_couple_couple_pr + Einter_couple_dipole_im &
            + Einter_couple_couple_im + Einter_couple_strain + Eel_strain
    E0 = Etotal + K0*log(rc) + K2/rc**2
    
    IF (verbosity.GE.verbosity_max) THEN
            IF (strain) THEN
                    WRITE(out,'(2(a,g16.8),a)') '  Homogeneous strain contribution: &
                        &                    E0 = ', Eel_strain, " = ", Eel_strain/length, " per unit length"
                    IF (l_dislo) WRITE(out,'(2(a,g16.8),a)') &
                        '                 interaction with dislocation dipoles: E0 = ', &
                        Einter_dipole_strain, ' = ', Einter_dipole_strain/length, ' per unit length'
                    IF (l_lineCouple) WRITE(out,'(2(a,g16.8),a)') &
                        '                 interaction with line-force couples:  E0 = ', &
                        Einter_couple_strain, ' = ', Einter_couple_strain/length, ' per unit length'
                    WRITE(out,*)
            END IF
            IF ((l_dislo.AND.(nd.GT.2)).OR.(l_lineCouple.AND.(nlc.GT.1)).OR.(l_dislo.AND.l_lineCouple)) THEN
                    WRITE(out,'(a)') '  Interaction between line-defects contained in primitive unit cell'
                    IF (l_dislo.AND.(nd.GT.2)) WRITE(out,'(2(a,g16.8),a)') &
                       '         dislocation dipoles - dislocation dipoles:    E0 = ', &
                       Einter_dipole_dipole_pr, ' = ', Einter_dipole_dipole_pr/length, ' per unit length'
                    IF (l_lineCouple.AND.(nlc.GT.1)) WRITE(out,'(2(a,g16.8),a)') &
                       '           line-force couples - line-force couples:    E0 = ', &
                       Einter_couple_couple_pr, ' = ', Einter_couple_couple_pr/length, ' per unit length'
                    IF (l_dislo.AND.l_lineCouple)  WRITE(out,'(2(a,g16.8),a)') &
                       '          dislocation dipoles - line-force couples:    E0 = ', &
                       Einter_dipole_couple_pr, ' = ', Einter_dipole_couple_pr/length, ' per unit length'
                    WRITE(out,*)
            END IF
            IF (l_dislo.AND.(xImages.OR.yImages.OR.zImages)) THEN
                    WRITE(out,'(a)') '  Dislocation interaction with image line-defects'
                    WRITE(out,'(2(a,g16.8),a)') '                                image dislocations:    E0 = ', &
                       Einter_dipole_dipole_im, ' = ', Einter_dipole_dipole_im/length, ' per unit length'
                    IF (l_lineCouple)  WRITE(out,'(2(a,g16.8),a)') &
                       '                          image line-force couples:    E0 = ', &
                       0.5d0*Einter_couple_dipole_im, ' = ', 0.5d0*Einter_couple_dipole_im/length, ' per unit length'
                    WRITE(out,*)
            END IF
            IF (l_lineCouple.AND.(xImages.OR.yImages.OR.zImages)) THEN
                    WRITE(out,'(a)') '  Line-force couple interaction with image line-defects'
                    WRITE(out,'(2(a,g16.8),a)') , '                          image line-force couples:    E0 = ', &
                       Einter_couple_couple_im, ' = ', Einter_couple_couple_im/length, ' per unit length'
                    IF (l_dislo)  WRITE(out,'(2(a,g16.8),a)') &
                       '                                image dislocations:    E0 = ', &
                       0.5d0*Einter_couple_dipole_im, ' = ', 0.5d0*Einter_couple_dipole_im/length, ' per unit length'
                    WRITE(out,*)
            END IF
            WRITE(out,'(a)') '--------------'
            WRITE(out,*)
            WRITE(out,'(2(a,g16.8),a)') ' Total elastic energy :    E = E0 - K0 ln(rc) - K2 / rc^2 = ', Etotal, ' = ', Etotal/length, ' per unit length'
            WRITE(out,'(2(a,g16.8),a)') '       with E0 = ', E0, ' = ', E0/length, ' per unit length'
            WRITE(out,'(2(a,g16.8),a)') '            K0 = ', K0, ' = ', K0/length, ' per unit length'
            WRITE(out,'(2(a,g16.8),a)') '            K2 = ', K2, ' = ', K2/length, ' per unit length'
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
    END IF

  END SUBROUTINE ElasticEnergy

END MODULE EnergyModule
