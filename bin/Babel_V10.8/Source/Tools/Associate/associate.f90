MODULE Associate_module

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE Associate_Atoms(out)

    USE MathStruct
    USE Associate_data
    IMPLICIT NONE

    INTEGER, intent(in) :: out

    INTEGER :: i, i0, i0min, nMove
    REAL(kind(0.d0)), dimension(1:3) :: x, s, s0, dr
    REAL(kind(0.d0)) :: dr2, dr2_min
    REAL(kind(0.d0)), dimension(:,:), allocatable :: xc, xc0
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at, inv_at0

    ! Initialization
    ind(:)=0 ; ind0(:)=0
    nMove = 0

    ! Check if periodic boundary conditions are used
    IF (periodic) THEN
            IF (Sum(at(:,:)**2).LE.distance_zero2) THEN
                     WRITE(0,'(a)') 'Basis vectors at(1:3,1:3) undefined for input structure'
                     WRITE(0,'(a)') '  => could not apply periodic boundary conditions'
                     WRITE(0,'(a)') '     You need to define at(1:3,1:3) or to set periodic to .false. '
                     STOP '< Associate_Atoms >'
            END IF
            IF (Sum(at0(:,:)**2).LE.distance_zero2) THEN
                     WRITE(0,'(a)') 'Basis vectors at0(1:3,1:3) undefined for reference structure'
                     WRITE(0,'(a)') '  => could not apply periodic boundary conditions'
                     WRITE(0,'(a)') '     You need to define at0(1:3,1:3) or to set periodic to .false. '
                     STOP '< Associate_Atoms >'
            END IF

            CALL Mat3Inv(at,inv_at)
            IF (Allocated(xc)) Deallocate(xc)
            ALLOCATE(xc(1:3,1:imm))
            xc(:,1:im)=MatMul(inv_at(:,:),xp(:,1:im))

            CALL Mat3Inv(at0,inv_at0)
            IF (Allocated(xc0)) Deallocate(xc0)
            ALLOCATE(xc0(1:3,1:imm))
            xc0(:,1:im0)=MatMul(inv_at0(:,:),xp0(:,1:im))

    END IF


    DO i=1, im    ! Loop on atoms in input structure
       x(:) = xp(:,i)
       IF (periodic) s(:) = xc(:,i)
       i0min=0
       dr2_min=-1.d0
       DO i0=1, im0    ! Loop on atoms in input structure
          IF (periodic) THEN
                  s0(:) = xc0(:,i0)
                  s = s - aNInt( s-s0 )
                  dr = MatMul(at,s) - MatMul(at0,s0) 
          ELSE
                  dr(:) =  x(:) - xp0(:,i0) 
          END IF
          dr2 = Sum( dr(:)**2 )
          IF ( (dr2_min.LT.0.d0).OR.(dr2.LT.dr2_min) ) THEN
                  dr2_min=dr2
                  i0min=i0
          END IF
       END DO
       IF ( (i0min.EQ.0).OR.(i0min.GT.im0)) THEN
               WRITE(0,'(a,i0)') 'Does not manage to find closest replica for atom ', i
               STOP '< Associate_Atoms >'
       ELSE IF (ind0(i0min).NE.0) THEN
               WRITE(0,'(2(a,i0))') 'Two atoms share the same replica: atoms ', i, &
                  ' and ', ind0(i0min)
               STOP '< Associate_Atoms >'
       ELSE
               ind(i)=i0min
               ind0(i0min)=i
               IF (i.NE.i0min) THEN
                       nMove = nMove + 1
                       IF (verbosity.LE.verbosity_max) THEN
                               WRITE(out,'(2(a,i0),a,g20.6)') &
                                        ' atom ', i, ' associated with atom ', i0min, &
                                        ':  dr^2=', dr2_min
                       END IF
               END IF
       END IF
    END DO

    IF (periodic) THEN
            DEALLOCATE(xc)
            DEALLOCATE(xc0)
    END IF

    IF (verbosity.LE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(i0,a)') nMove, ' atoms have been displaced'
            WRITE(out,*)
    END IF


  END SUBROUTINE Associate_Atoms

  ! ==============================================

  SUBROUTINE Sort_Atoms(xp, iTyp, im, ind)
 
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(:,:), intent(inout) :: xp
    INTEGER, dimension(:), intent(inout) :: iTyp
    INTEGER, dimension(:), intent(inout) :: ind
    INTEGER, intent(in) :: im

    REAL(kind(0.d0)), dimension(1:3,1:im) :: xpn
    INTEGER, dimension(1:im) :: iTypn
    INTEGER :: i, in

    DO i=1, im
      in=ind(i)
      xpn(:,in) = xp(:,i)
      iTypn(in) = iTyp(i)
      ind(i)=i
    END DO

    xp(1:3,1:im) = xpn(1:3,1:im)
    iTyp(1:im) = iTypn(1:im)
  

  END SUBROUTINE Sort_Atoms

END MODULE Associate_module
