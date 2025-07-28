MODULE PeriodRelaxModule

  USE Babel_data

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE RelaxPeriod(epsi, out)

    USE lineCouple_elasticity_ani, ONLY : l_lineCouple
    USE LineForce_elasticity_ani, ONLY : lineForce
    USE disloc_elasticity_ani, ONLY : l_dislo
    USE DDipoleModule, ONLY : l_DDipole
    USE loopModule, ONLY : l_loop
    IMPLICIT NONE
    INTEGER, intent(in), optional :: out  ! output unit
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: epsi

    REAL(kind(0.d0)), dimension(1:3,1:3) :: epsi1

    epsi(:,:) = 0.d0

    IF (l_dislo) THEN
            IF (Present(out)) THEN
                    CALL RelaxPeriod_dislo(epsi1,out)
            ELSE
                    CALL RelaxPeriod_dislo(epsi1)
            ENDIF
            epsi = epsi + epsi1
    END IF
    IF (l_DDipole) THEN
            IF (Present(out)) THEN
                    CALL RelaxPeriod_DDipole(epsi1,out)
            ELSE
                    CALL RelaxPeriod_DDipole(epsi1)
            ENDIF
            epsi = epsi + epsi1
    END IF
    IF (l_loop) THEN
            IF (Present(out)) THEN
                    CALL RelaxPeriod_Loops(epsi1,out)
            ELSE
                    CALL RelaxPeriod_Loops(epsi1)
            ENDIF
            epsi = epsi + epsi1
    END IF
    IF (lineForce) THEN
            IF (Present(out)) THEN
                    CALL RelaxPeriod_lineForce(epsi1,out)
            ELSE
                    CALL RelaxPeriod_lineForce(epsi1)
            ENDIF
            epsi = epsi + epsi1
    END IF
    IF (l_lineCouple) THEN
            IF (Present(out)) THEN
                    CALL RelaxPeriod_lineCouple(epsi1,out)
            ELSE
                    CALL RelaxPeriod_lineCouple(epsi1)
            ENDIF
            epsi = epsi +  epsi1
    END IF

  END SUBROUTINE RelaxPeriod

  !====================================================================

  SUBROUTINE RelaxPeriod_dislo(epsi,out)

    USE babel_data
    USE disloc_elasticity_ani
    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: epsi
    INTEGER, intent(in), optional :: out  ! output unit

    REAL(kind(0.d0)), dimension(1:3) :: area
    REAL(kind(0.d0)) :: volume, height, inv_area
    INTEGER :: i, j, n

    ! Initialization
    epsi=0.d0

    IF (xImages.OR.yImages.OR.zImages) THEN ! Periodic boundary conditions

            ! Displacement gradient due to dislocation distribution
            epsi = 0.d0
            DO n=1, nd
               area(1:3) = CrossProduct( lDislo(1:3,n), cDislo(1:3,n) )
               DO j=1,3
                  epsi(1:3,j) = epsi(1:3,j) + bDislo(1:3,n)*area(j)
               END DO
            END DO

            ! Normalize displacement gradient by unit cell area
            Volume = Abs( ScalarTripleProduct( at(1:3,1), at(1:3,2), at(1:3,3)))
            IF (yImages.AND.zImages) THEN
                    height = Sqrt( Sum( at(1:3,1)**2 ) )
            ELSE IF (xImages.AND.zImages) THEN
                    height = Sqrt( Sum( at(1:3,2)**2 ) )
            ELSE IF (xImages.AND.yImages) THEN
                    height = Sqrt( Sum( at(1:3,3)**2 ) )
            ELSE
                    ! We have infinite dimension in one of the direction
                    ! orthogonal to the line direction of the dislocation
                    ! distribution
                    height = 0.d0
            END IF
            inv_area = height/volume
            epsi = inv_area*epsi

            IF (verbosity.GE.verbosity_max) THEN
                    IF (.NOT.Present(out)) Return
                    WRITE(out,*)
                    WRITE(out,'(a)') ' Displacement gradient equivalent to&
                        & dislocation distribution'
                    DO i=1,3
                       WRITE(out,'(a,i1,a,3(g11.4,2x),a)') '  e(',i,',1:3) = | ', &
                            epsi(i,1:3), ' |'
                    END DO
                    WRITE(out,'(a,g22.14)') '   area = ', 1.d0/inv_area
                    WRITE(out,*)
            END IF

    END IF

  END SUBROUTINE RelaxPeriod_dislo

  !====================================================================

  SUBROUTINE RelaxPeriod_DDipole(epsi,out)

    USE babel_data
    USE DDipoleModule
    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: epsi
    INTEGER, intent(in), optional :: out  ! output unit

    REAL(kind(0.d0)), dimension(1:3) :: area
    REAL(kind(0.d0)) :: volume, height, inv_area
    INTEGER :: i, j, n

    ! Initialization
    epsi=0.d0

    IF (xImages.OR.yImages.OR.zImages) THEN ! Periodic boundary conditions

            ! Displacement gradient due to dislocation distribution
            epsi = 0.d0
            DO n=1, nDDipole
               area(1:3) = CrossProduct( lDDipole(1:3,n), c1DDipole(1:3,n)-c2DDipole(1:3,n) )
               DO j=1,3
                  epsi(1:3,j) = epsi(1:3,j) + bDDipole(1:3,n)*area(j)
               END DO
            END DO

            ! Normalize displacement gradient by unit cell area
            Volume = Abs( ScalarTripleProduct( at(1:3,1), at(1:3,2), at(1:3,3)))
            IF (yImages.AND.zImages) THEN
                    height = Sqrt( Sum( at(1:3,1)**2 ) )
            ELSE IF (xImages.AND.zImages) THEN
                    height = Sqrt( Sum( at(1:3,2)**2 ) )
            ELSE IF (xImages.AND.yImages) THEN
                    height = Sqrt( Sum( at(1:3,3)**2 ) )
            ELSE
                    ! We have infinite dimension in one of the direction
                    ! orthogonal to the line direction of the dislocation
                    ! distribution
                    height = 0.d0
            END IF
            inv_area = height/volume
            epsi = inv_area*epsi

            IF (verbosity.GE.verbosity_max) THEN
                    IF (.NOT.Present(out)) Return
                    WRITE(out,*)
                    WRITE(out,'(a)') ' Displacement gradient equivalent to&
                        & dislocation dipole distribution'
                    DO i=1,3
                       WRITE(out,'(a,i1,a,3(g11.4,2x),a)') '  e(',i,',1:3) = | ', &
                            epsi(i,1:3), ' |'
                    END DO
                    WRITE(out,'(a,g22.14)') '   area = ', 1.d0/inv_area
                    WRITE(out,*)
            END IF

    END IF

  END SUBROUTINE RelaxPeriod_DDipole

  !====================================================================

  SUBROUTINE RelaxPeriod_Loops(epsi,out)

    USE babel_data
    USE loopModule
    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: epsi
    INTEGER, intent(in), optional :: out  ! output unit

    REAL(kind(0.d0)), dimension(1:3) :: area, xc, rA, rB
    REAL(kind(0.d0)) :: volume
    INTEGER :: i, il, j, nl

    ! Initialization
    epsi=0.d0


            ! Displacement gradient due to dislocation loops
            DO nl=1, nLoop

               ! Loop center 
               xc(:) = cLoop(:,nl)

               area(:) = 0.d0
               DO il=1, iLoop(nl)
                  ! Add area of triangular loop ABC
                  rA(:) = xLoop(:,il-1,nl) - xC(:)
                  rB(:) = xLoop(:,il,nl) - xC(:)
                  area(:) = area(:) - 0.5d0*CrossProduct(rA,rB)
               END DO

               DO j=1,3
                  epsi(1:3,j) = epsi(1:3,j) + bLoop(1:3,nl)*area(j)
               END DO
               
            END DO

            ! Normalize displacement gradient by unit cell area
            Volume = Abs( ScalarTripleProduct( at(1:3,1), at(1:3,2), at(1:3,3)))
            epsi(:,:) = epsi(:,:)/volume

            IF (verbosity.GE.verbosity_max) THEN
                    IF (.NOT.Present(out)) Return
                    WRITE(out,*)
                    WRITE(out,'(a)') ' Displacement gradient equivalent to&
                        & loop population'
                    DO i=1,3
                       WRITE(out,'(a,i1,a,3(g11.4,2x),a)') '  e(',i,',1:3) = | ', &
                            epsi(i,1:3), ' |'
                    END DO
                    WRITE(out,*)
            END IF


  END SUBROUTINE RelaxPeriod_Loops

  !====================================================================

  SUBROUTINE RelaxPeriod_lineForce(epsi,out)

    USE babel_data
    USE LineForce_elasticity_ani
    USE elasticity_ani
    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: epsi
    INTEGER, intent(in), optional :: out  ! output unit

    ! Elastic constants in initial and line-force axes
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: elastic_C, elastic_Cp
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CpVoigt, SpVoigt
    REAL(kind(0.d0)), dimension(1:6) :: sVoigt, eVoigt

    REAL(kind(0.d0)), dimension(1:3,1:3) :: forceDipole, rotDistrib, inv_rotDistrib
    REAL(kind(0.d0)), dimension(1:3) :: ex, ey, ez
    REAL(kind(0.d0)) :: ex_norm2, inv_norm, inv_area
    INTEGER :: i, j, n

    ! Initialization
    epsi=0.d0

    IF (xImages.OR.yImages.OR.zImages) THEN ! Periodic boundary conditions

            elastic_C=Unpack_CVoigt(CVoigt)

            ! Calculate dipolar tensor equivalent to line force distribution
            forceDipole(1:3,1:3) = 0.d0
            DO n=1, nlf
               DO j=1,3
                  forceDipole(1:3,j) = forceDipole(1:3,j) + fLineForce(1:3,n)*cLineForce(j,n)
               END DO
            END DO

            ! Orientate the crystal
            ez(1:3) = lLineForce(1:3,1)      ! Line direction of the distribution
            ex(1)=0.d0 ; ex(2)=ez(3) ; ex(3)=-ez(2)     ! ez ^ (1,0,0)
            ex_norm2 = Sum( ex(1:3)**2 )
            IF (ex_norm2.LE.distance_Zero2) THEN
                    ex(1)=-ez(3) ; ex(2)=0.d0 ; ex(3)=ez(1)     ! ez ^ (0,1,0)
                    ex_norm2 = Sum( ex(1:3)**2 )
            END IF
            inv_norm = 1.d0/Sqrt( ex_norm2 )
            ex(1:3) = inv_norm*ex(1:3)
            ey(1:3) = CrossProduct( ez(1:3), ex(1:3) )
            rotDistrib(1,1:3) = ex(1:3)
            rotDistrib(2,1:3) = ey(1:3)
            rotDistrib(3,1:3) = ez(1:3)
            inv_rotDistrib = Transpose(rotDistrib)

            ! Rotate elastic constants in the corresponding axes
            elastic_Cp = Rotate_Elastic_Constants(elastic_C, rotDistrib(1:3,1:3) )
            CpVoigt = Pack_CVoigt(elastic_Cp)
            CALL MatInv(CpVoigt, SpVoigt)

            ! Rotate dipolar tensor
            forceDipole = MatMul( MatMul( rotDistrib, forceDipole ), inv_rotDistrib )

            ! Equivalent stress
            sVoigt(1:6) = (/ forceDipole(1,1), forceDipole(2,2), forceDipole(3,3), &
                0.5d0*(forceDipole(2,3)+forceDipole(3,2)), &
                0.5d0*(forceDipole(1,3)+forceDipole(3,1)), &
                0.5d0*(forceDipole(1,2)+forceDipole(2,1)) /)

            ! Displacement gradient in dipolar axes
            eVoigt = MatMul(SpVoigt, sVoigt)
            epsi = ReShape( (/ eVoigt(1), 0.5d0*eVoigt(6), 0.5d0*eVoigt(5), &
                                      0.5d0*eVoigt(6),       eVoigt(2), 0.5d0*eVoigt(4), &
                                      0.5d0*eVoigt(5), 0.5d0*eVoigt(4),       eVoigt(3) /), &
                                      (/3,3/) )

            ! Rotate displacement gradient
            epsi = MatMul( MatMul( inv_rotDistrib, epsi ), rotDistrib )
            
            ! Normalize displacement gradient by unit cell area
            IF (yImages.AND.zImages) THEN
                    ex = CrossProduct( at(1:3,2), at(1:3,3) ) 
                    inv_area = 1.d0/Sqrt( Sum( ex(1:3)**2 ) )
            ELSE IF (xImages.AND.zImages) THEN
                    ey = CrossProduct( at(1:3,1), at(1:3,3) )
                    inv_area = 1.d0/Sqrt( Sum( ey(1:3)**2 ) )
            ELSE IF (xImages.AND.yImages) THEN
                    ez = CrossProduct( at(1:3,1), at(1:3,2) )
                    inv_area = 1.d0/Sqrt( Sum( ez(1:3)**2 ) )
            ELSE
                    ! We have infinite dimension in one of the direction
                    ! orthogonal to the line direction of the dislocation
                    ! distribution
                    inv_area = 0.d0
            END IF
            epsi = inv_area*epsi

            IF (verbosity.GE.verbosity_max) THEN
                    IF (.NOT.Present(out)) Return
                    WRITE(out,*)
                    WRITE(out,'(a)') ' Dipolar tensor equivalent to line force distribution '
                    DO i=1,3
                       WRITE(out,'(a,3g14.6,a)') '   | ', forceDipole(i,1:3), ' |'
                    END DO
                    WRITE(out,'(a,3g14.6)') '  with line direction: ', ez(1:3)
                    WRITE(out,*)
                    WRITE(out,'(a)') '  displacement gradient equivalent to&
                        & line force distribution'
                    DO i=1,3
                       WRITE(out,'(a,i1,a,3(g11.4,2x),a)') '  e(',i,',1:3) = | ', &
                            epsi(i,1:3), ' |'
                    END DO
                    WRITE(out,*)
            END IF

    END IF

  END SUBROUTINE RelaxPeriod_lineForce

  !====================================================================

  SUBROUTINE RelaxPeriod_lineCouple(epsi,out)

    USE babel_data
    USE elasticity_ani
    USE elasticity_Stroh
    USE LineCouple_elasticity_ani
    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(out) :: epsi
    INTEGER, intent(in), optional :: out  ! output unit

    ! Elastic constants in initial and line-force axes
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: elastic_C, elastic_Cp

    REAL(kind(0.d0)), dimension(1:3,1:3) :: epsi1
    REAL(kind(0.d0)) :: inv_area
    REAL(kind(0.d0)), dimension(1:3) :: ex, ey, ez
    INTEGER :: n, i
    REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
    REAL(kind(0.d0)), dimension(1:6) :: sVoigt, eVoigt
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CpVoigt, SpVoigt


    epsi(:,:) = 0.d0
            
    IF (xImages.OR.yImages.OR.zImages) THEN ! Periodic boundary conditions

            elastic_C=Unpack_CVoigt(CVoigt)
            
            loop_couple: DO n=1, nlc

            ! Rotate elastic constants in the line-force axes
            elastic_Cp = Rotate_Elastic_Constants(elastic_C, rotLineCouple(1:3,1:3,n) )
            
                    ! Calculate unit cell area  for displacement gradient
                    IF (yImages.AND.zImages) THEN
                            ex = CrossProduct( at(1:3,2), at(1:3,3) ) 
                            inv_area = 1.d0/Sqrt( Sum( ex(1:3)**2 ) )
                    ELSE IF (xImages.AND.zImages) THEN
                            ey = CrossProduct( at(1:3,1), at(1:3,3) )
                            inv_area = 1.d0/Sqrt( Sum( ey(1:3)**2 ) )
                    ELSE IF (xImages.AND.yImages) THEN
                    ez = CrossProduct( at(1:3,1), at(1:3,2) )
                    inv_area = 1.d0/Sqrt( Sum( ez(1:3)**2 ) )
                    ELSE
                            ! We have infinite dimension in one of the direction
                            ! orthogonal to the line direction 
                            inv_area = 0.d0
                    END IF

                    ! Homogeneous stress
                    !sVoigt(1:6) = (/ mxLineCouple(n)*inv_area, myLineCouple(n)*inv_area, 0.d0, &
                    !    0.d0, 0.d0, 0.d0 /)
                    sVoigt(1:6) = MatMul( inter_lineCouple_Voigt(1:6,1:2, n), &
                         (/ mxLineCouple(n)*inv_area, myLineCouple(n)*inv_area /) )
                    ! Voigt notation of rotated elastic constants
                    CpVoigt = Pack_CVoigt(elastic_Cp)
                    ! Inverse matrix (compliances)
                    CALL MatInv(CpVoigt, SpVoigt)
                    ! Homogeneous strain
                    eVoigt(:) = MatMul( SpVoigt, sVoigt )
                    ! Displacement in cartesian coordinates
                    epsi1(:,:) = MatMul( inv_rotLineCouple(1:3,1:3,n), &
                        MatMul( ReShape( (/ eVoigt(1), 0.5d0*eVoigt(6), 0.5d0*eVoigt(5), &
                                      0.5d0*eVoigt(6),       eVoigt(2), 0.5d0*eVoigt(4), &
                                      0.5d0*eVoigt(5), 0.5d0*eVoigt(4),       eVoigt(3) /), &
                                      (/3,3/) ), rotLineCouple(1:3,1:3,n) ) )
                     ! Total strain
                     epsi(:,:) = epsi(:,:) + epsi1(:,:)
                    
                   IF (verbosity.GE.verbosity_max) THEN
                           IF (.NOT.Present(out)) Cycle
                           IF (inv_area.NE.0.d0) THEN
                                   WRITE(out,'(a,i0)') '  displacement gradient equivalent to line-force couple ', n
                                   DO i=1,3
                                      WRITE(out,'(a,i1,a,3g14.6)') '  e(',i,',1:3) = ', epsi1(i,1:3)
                                   END DO
                           END IF
                           WRITE(out,*)
                   END IF

            END DO loop_couple
    END IF
    
  END SUBROUTINE RelaxPeriod_lineCouple

END MODULE PeriodRelaxModule
