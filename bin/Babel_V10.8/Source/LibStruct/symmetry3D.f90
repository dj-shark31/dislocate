MODULE Symmetry3DModule

        IMPLICIT NONE
        
        ! Maximal number of symmetry operations
        !   (without translations corresponding to PBC)
        INTEGER, parameter :: max_nSym=48

        ! Symmetry operation
        TYPE sym3D_op_t
                ! Rotation part in crystal and cartesian coordinates
                INTEGER, dimension(1:3,1:3) :: iRot
                REAL(kind(0.d0)), dimension(1:3,1:3) :: rot
                ! Translation part in crystal and cartesian coordinates
                REAL(kind(0.d0)), dimension(1:3) :: uReduced
                REAL(kind(0.d0)), dimension(1:3) :: u
        END TYPE sym3D_op_t

        ! Symmetry group
        TYPE sym3D_group_t
                ! Number of symmetry operations
                INTEGER :: nSym 
                ! Symmetry operations
                TYPE(sym3D_op_t), dimension(1:max_nSym) :: op
                ! Translations due to periodic boundary condtions
                !   periodicity vectors:  uPBC(1:3,1), uPBC(1:3,2), uPBC(1:3,3)
                REAL(kind(0.d0)), dimension(1:3,1:3) :: uPBC
        END TYPE sym3D_group_t

        ! Zero for distances (in A)
        REAL(kind(0.d0)), parameter, private :: Zero=1.d-4
        REAL(kind(0.d0)), parameter, private :: Zero2=Zero*Zero

CONTAINS

        SUBROUTINE InitSym3D(sym)
                ! Initialize symmetry operation sym

                IMPLICIT NONE
                TYPE(sym3D_op_t), intent(out) :: sym

                sym%rot(:,:) = 0.d0  ; sym%iRot(:,:) = 0
                sym%rot(1,1) = 1.d0  ; sym%iRot(1,1) = 1
                sym%rot(2,2) = 1.d0  ; sym%iRot(2,2) = 1
                sym%rot(3,3) = 1.d0  ; sym%iRot(3,3) = 1
                sym%u(:) = 0.d0      ; sym%uReduced(:) = 0.d0
        
        END SUBROUTINE InitSym3D

        SUBROUTINE InitSym3Dgroup(group)
                ! Initialize symmetry operation group group

                IMPLICIT NONE 
                TYPE(sym3D_group_t), intent(out) :: group

                INTEGER :: n

                group%nSym=0
                group%uPBC(:,:) = 0.d0
                DO n=1, max_nSym
                   CALL InitSym3D(group%op(n))
                END DO

        END SUBROUTINE InitSym3Dgroup

        FUNCTION CheckBasisSym3D(at,inv_at, sym) RESULT(test)
                ! Check that periodicity vectors at 
                !   are compatible with symmetry operation sym
  
                IMPLICIT NONE
                ! Periodicity vectors: at(1:3,i) and inverse
                REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at, inv_at
                ! Symmetry operation
                TYPE(sym3D_op_t), intent(in) :: sym
  
                INTEGER :: i
                REAL(kind(0.d0)), dimension(1:3) :: atNew, s, ds
                LOGICAL :: test
  
                test=.true.
                DO i=1, 3
                      ! Apply symmetry operations
                      atNew(1:3) = MatMul(sym%rot(1:3,1:3), at(1:3,i) )
                      ! Reduced units
                      s(1:3) = MatMul( inv_at(1:3,1:3), atNew(1:3) )
                      ds(1:3) = s(1:3) - aNInt(s(1:3))
                      IF ( (Abs(ds(1)).GT.zero).OR.(Abs(ds(2)).GT.zero).OR.(Abs(ds(3)).GT.zero) ) &
                              test=.FALSE.
                END DO
  
        END FUNCTION CheckBasisSym3D

        SUBROUTINE CheckBasisSym3Dgroup(at, group, out)
                 
                USE MathStruct
                IMPLICIT NONE
                ! Periodicity vectors: at(1:3,i)
                REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
                ! Group of symmetry operations
                TYPE(sym3D_group_t), intent(in) :: group
                ! Output unit
                INTEGER, intent(in) :: out

                REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
                LOGICAL :: test 
                INTEGER :: ns

                ! Periodicity vectors
                CALL Mat3inv(at, inv_at)
                test=.TRUE.
                DO ns=1, group%nSym
                   IF (.NOT.CheckBasisSym3D(at, inv_at, group%op(ns))) THEN
                           WRITE(out,'(a,i0)') 'Periodicity vectors&
                                      & are not compatible with symetry operation ', ns
                           CALL PrintSym3D(group%op(ns),out)
                           test=.FALSE.
                   END IF
                END DO
                IF (.NOT.test) THEN
                        WRITE(out,'(a)') '  Periodicity vectors:'
                        WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
                        WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
                        WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
                        STOP '< CheckBasisSym3Dgroup >'
                END IF

        END SUBROUTINE CheckBasisSym3Dgroup

        SUBROUTINE CheckSym3Dgroup(group, out)
                 
                IMPLICIT NONE
                ! Group of symmetry operations
                TYPE(sym3D_group_t), intent(in) :: group
                ! Output unit
                INTEGER, intent(in) :: out

                ! Check periodicity vectors defined in symmetry group
                CALL CheckBasisSym3Dgroup(group%uPBC, group, out)

        END SUBROUTINE CheckSym3Dgroup

        FUNCTION EquivalentSymmetry3D(sym1, sym2) RESULT(test)
                ! Function returns .true. or .false. if the symmetry operations are
                ! equivalent or not

                IMPLICIT NONE
                TYPE(sym3D_op_t), intent(in) :: sym1, sym2
                LOGICAL :: test

                INTEGER :: irotNorm2
                REAL(kind(0.d0)) :: uNorm2
                REAL(kind(0.d0)), dimension(1:3) :: dur

                ! Compare rotation matrices
                irotNorm2 = Sum( ( sym1%irot(:,:) - sym2%irot(:,:) )**2 )
                IF (irotNorm2.GT.0) THEN
                        test=.FALSE.
                        RETURN
                END IF

                ! Compare translation vectors, taking into account periodicity
                dur(:) = sym1%uReduced(:) - sym2%uReduced(:)
                dur(:) = dur(:) - aNInt( dur(:) )
                uNorm2 = Sum( dur(:)**2 )
                IF (uNorm2.GT.zero2) THEN
                        test=.FALSE.
                        RETURN
                END IF

                test=.TRUE.
                RETURN

        END FUNCTION EquivalentSymmetry3D

        FUNCTION InsideSpaceGroup3D(sym0, symGroup) RESULT(test)
                ! Returns .true. if symmetry operation sym0 belongs to sym2D(1:nSym)

                IMPLICIT NONE
                TYPE(sym3D_op_t) :: sym0
                TYPE(sym3D_group_t), intent(in) :: symGroup
                LOGICAL :: test

                INTEGER :: ns

                DO ns=1, symGroup%nSym
                   IF ( EquivalentSymmetry3D(sym0, symGroup%op(ns)) ) THEN
                           test=.TRUE.
                           RETURN
                   END IF
                END DO
                test=.FALSE.
                RETURN

        END FUNCTION InsideSpaceGroup3D

        FUNCTION ComposeSymmetry3D( sym1, sym2) RESULT(sym3)
                ! Generate symmetry operation sym3 by composing sym1 and sym2

                IMPLICIT NONE
                ! Input symmetry operations
                TYPE(sym3D_op_t), intent(in) :: sym1, sym2
                ! Output symmetry operation
                TYPE(sym3D_op_t) :: sym3

                sym3%rot(:,:) = MatMul( sym1%rot(:,:), sym2%rot(:,:) )
                sym3%irot(:,:) = MatMul( sym1%irot(:,:), sym2%irot(:,:) )
                sym3%u(:) = sym1%u(:) + MatMul( sym1%rot(:,:), sym2%u(:) )
                sym3%uReduced(:) = sym1%uReduced(:) + MatMul( Dble( sym1%irot(:,:) ), sym2%uReduced(:) )

        END FUNCTION ComposeSymmetry3D

        SUBROUTINE BuildSpaceGroup3D(symGroup0, symGroup)
                ! Builds space group symGroup containing the input
                ! symmetry operations symGroup0

                IMPLICIT NONE
                ! Input symmetry operations
                TYPE(sym3D_group_t), intent(in) :: symGroup0
                ! Output symmetry operations (space group)
                TYPE(sym3D_group_t), intent(out) :: symGroup

                INTEGER :: ns0, ns1, ns2
                TYPE(sym3D_op_t) :: sym3

                ! Periodicity vectors
                symGroup%uPBC = symGroup0%uPBC

                ! First symmetry operations in point group = identity
                symGroup%nSym = 1
                symGroup%op(1)%rot(:,:)  = Reshape( (/ 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0 /), (/3,3/) )
                symGroup%op(1)%irot(:,:) = Reshape( (/ 1,0,0, 0,1,0, 0,0,1 /), (/3,3/) )
                symGroup%op(1)%u(:) = 0.d0
                symGroup%op(1)%uReduced(:) = 0.d0

                ! Add one instance of each input symmetry operation
                DO ns0=1, symGroup0%nSym
                   IF (.NOT.InsideSpaceGroup3D(symGroup0%op(ns0), symGroup)) THEN
                           symGroup%nSym = symGroup%nSym + 1
                           symGroup%op(symGroup%nSym) = symGroup0%op(ns0)
                           !WRITE(6,'(a,i0,a)') 'Structure ', ns0, ' added'          ! DEBUG
                   !ELSE                                                             ! DEBUG
                           !WRITE(6,'(a,i0,a)') 'Structure ', ns0, ' already inside' ! DEBUG
                   END IF
                END DO

                ! Generate new symmetry operations by composing them and add them to the
                ! space group if they are not already included
                ns1 = 2
                loop_sym1: DO

                   ns2 = 2
                   loop_sym2: DO
                    
                      sym3 =  ComposeSymmetry3D( symGroup%op(ns1), symGroup%op(ns2))  
                      !WRITE(6,*)                                                               ! DEBUG
                      !WRITE(6,'(2(a,i0))') "Result of symmetry operations ", ns1, " and ", ns2 ! DEBUG
                      !CALL PrintSym3D(sym3, 6)                                            ! DEBUG
                      IF (.NOT.InsideSpaceGroup3D(sym3, symGroup)) THEN
                              symGroup%nSym = symGroup%nSym + 1
                              symGroup%op(symGroup%nSym) = sym3
                              !WRITE(6,'(a,i0)') 'Structure added as number ', symGroup%nSym ! DEBUG
                      !ELSE                                                         ! DEBUG
                              !WRITE(6,'(a)') 'Structure already inside'            ! DEBUG
                      END IF

                      ! Increment symmetry operation 
                      ns2 = ns2 + 1
                      IF (ns2.GT.symGroup%nSym) EXIT
                   END DO loop_sym2

                   ! Increment symmetry operation 
                   ns1 = ns1 + 1
                   IF (ns1.GT.symGroup%nSym) EXIT

                END DO loop_sym1

        END SUBROUTINE BuildSpaceGroup3D  

        SUBROUTINE PrintSym3D(sym, out)
                ! Print on output unit out symmetry operation sym

                IMPLICIT NONE
                TYPE(sym3D_op_t), intent(in) :: sym
                INTEGER, intent(in) :: out

                WRITE(out,'(a,3g14.6,a,3i3,a)') '  rot(1,1:3) = (', sym%rot(1,1:3), ')   =  (', sym%irot(1,1:3), ') &
                        & (cart / crystal)'
                WRITE(out,'(a,3g14.6,a,3i3,a)') '  rot(2,1:3) = (', sym%rot(2,1:3), ')   =  (', sym%irot(2,1:3), ')'
                WRITE(out,'(a,3g14.6,a,3i3,a)') '  rot(3,1:3) = (', sym%rot(3,1:3), ')   =  (', sym%irot(3,1:3), ')'
                WRITE(out,'(a,3g14.6,a)')       '  u(1:3) = ', sym%u(1:3), ' (cartesian units)'
                WRITE(out,'(a,3g14.6,a)')       '         = ', sym%uReduced(1:3), ' (crystal units)'
                WRITE(out,*)

        END SUBROUTINE PrintSym3D

        SUBROUTINE PrintSym3Dgroup(group, out)
                ! Print on output unit out group of symmetry operations group

                IMPLICIT NONE
                TYPE(sym3D_group_t), intent(in) :: group
                INTEGER, intent(in) :: out

                INTEGER :: ns

                WRITE(out,*)
                WRITE(out,'(a)') ' --- Symmetry operations -----------------------'
                WRITE(out,*)
                WRITE(out,'(a,i0)') 'Number of symmetry operations, nSym = ', group%nSym
                DO ns=1, group%nSym
                   WRITE(out,*)
                   WRITE(out,'(a,i3)')  'Symmetry operation ', ns
                   CALL PrintSym3D(group%op(ns),out)
                END DO
                WRITE(out,*)
                WRITE(out,'(a)') 'Periodicity vectors'
                WRITE(out,'(2x,3g14.6)') group%uPBC(1:3,1)
                WRITE(out,'(2x,3g14.6)') group%uPBC(1:3,2)
                WRITE(out,'(2x,3g14.6)') group%uPBC(1:3,3)

                WRITE(out,*)
                WRITE(out,'(a)') ' --- End of Symmetry operations ----------------'
                WRITE(out,*)

        END SUBROUTINE PrintSym3Dgroup

        SUBROUTINE ApplySym3DReduced(xc0, xc1, im, sym)
                ! Apply symmetry operation sym to the atomic
                ! positions defined in crystal coordinates by xc0
                ! Result: xc1

                IMPLICIT NONE
                REAL(kind(0.d0)), dimension(:,:), intent(in)  :: xc0
                REAL(kind(0.d0)), dimension(:,:), intent(out) :: xc1
                INTEGER, intent(in) :: im
                TYPE(sym3D_op_t), intent(in) :: sym

                INTEGER :: n

                xc1(:,:) = 0.d0
                DO n=1, im
                   xc1(1:3,n) = MatMul( Dble( sym%irot(1:3,1:3) ), xc0(1:3,n) ) + sym%uReduced(1:3)
                END DO

        END SUBROUTINE ApplySym3DReduced
                
        SUBROUTINE ApplySym3D(xp0, xp1, im, sym)
                ! Apply symmetry operation sym to the atomic
                ! positions defined in cartesian coordinates by xp0
                ! Result: xp1

                IMPLICIT NONE
                REAL(kind(0.d0)), dimension(:,:), intent(in)  :: xp0
                REAL(kind(0.d0)), dimension(:,:), intent(out) :: xp1
                INTEGER, intent(in) :: im
                TYPE(sym3D_op_t), intent(in) :: sym

                INTEGER :: n

                xp1(:,:) = 0.d0
                DO n=1, im
                   xp1(1:3,n) = MatMul( sym%rot(1:3,1:3), xp0(1:3,n) ) + sym%u(1:3)
                END DO

        END SUBROUTINE ApplySym3D
                
        SUBROUTINE Symmetrize3DPositions(xp, iTyp, im, at, group, out, verbose)

                USE MathStruct
                IMPLICIT NONE
                ! Atom coordinates: xp(1:3,:)
                REAL(kind(0.d0)), dimension(:,:), intent(inout) ::xp
                ! Atom type
                INTEGER, dimension(:), intent(in) :: iTyp
                ! Number of atoms in simulation box
                INTEGER, intent(in) :: im
                ! Lattice vector coordinates (A): at(1:3,1), at(2:3,2), ...
                REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
                ! Symmetry operations
                TYPE(sym3D_group_t), intent(in) :: group
                ! Output unit
                INTEGER, intent(in) :: out
                ! Verbose mode
                LOGICAL, intent(in) :: verbose

                REAL(kind(0.d0)), dimension(1:3,1:im) :: xp0, xp1, xpS
                REAL(kind(0.d0)), dimension(1:3) :: dR, dS
                REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
                REAL(kind(0.d0)) :: at_norm2, dR2, dR2min
                INTEGER :: ns, i0, i1, iMin
                LOGICAL, dimension(1:im) :: taken
                LOGICAL :: at_defined

                IF (group%nSym.LE.0) THEN
                        WRITE(0,'(a)') 'You need to define symmetry operations'
                        WRITE(0,'(a,i0)') '  nSym = ', group%nSym
                        STOP '< SymmetrizePositions >'
                END IF

                ! Initialization
                xp0(1:3,1:im) = xp(1:3,1:im)    ! Initial positions
                xpS(:,:) = 0.d0                 ! Symmetrized positions

                ! Invert periodicity vectors if defined
                at_norm2=SUM( at(1:3,1:3)**2 )
                at_defined = (at_norm2.GT.zero2) 
                IF (at_defined) CALL Mat3Inv(at,inv_at)

                ! Loop on symmetry operations
                loop_sym: DO ns=1, group%nSym

                   ! Apply symmetry operation
                   CALL ApplySym3D(xp0, xp1, im, group%op(ns))
                
                   ! For each atom, look for closest atom of same type in initial structure
                   taken(:)=.FALSE.
                   loop_i1: DO i1=1, im

                      loop_i0: DO i0=1, im
                
                         IF (iTyp(i0).NE.iTyp(i1)) CYCLE
                
                         ! Vector linking both atoms
                         dR(:) = xp1(:,i1) - xp0(:,i0)

                         ! Apply periodic boundary conditions
                         IF (at_defined) THEN
                                 dS(:) = MatMul( inv_at(:,:), dR(:) )
                                 dS(:) = dS(:) - aNint( dS(:) )
                                 dR(:) = MatMul( at(:,:), dS(:) )
                         END IF

                         ! Check if the atom minimizes the distance
                         dR2 = Sum( dR(:)**2 )
                         IF ( (i0.EQ.1).OR.(dR2.LE.dR2min) ) THEN
                                 dR2min=dR2
                                 iMin=i0
                         END IF

                      END DO loop_i0

                      ! Check that atom has not already been chosen
                      IF (taken(iMin)) THEN
                              WRITE(0,'(a,i0)') 'Problem with symmetry ', ns
                              WRITE(0,'(a)')    '  program found two close atoms' 
                              STOP '< SymmetrizePositions >'
                      ELSE 
                              ! Add necessary shift corresponding to periodic boundary conditions
                              IF (at_defined) THEN
                                         dR(:) = xp1(:,i1) - xp0(:,iMin)
                                         dS(:) = MatMul( inv_at(:,:), dR(:) )
                                         dS(:) = aNint( dS(:) )
                                         dR(:) = MatMul( at(:,:), dS(:) )
                                         xp1(:,i1 ) = xp1(:,i1) - dR(:)
                              END IF

                              taken(iMin)=.TRUE.
                              xps(1:3,iMin) = xps(1:3,iMin) + xp1(1:3,i1)
                      END IF
                   END DO loop_i1

                END DO loop_sym

                ! Result
                xp(1:3,1:im) = xpS(1:3,1:im)/dble(group%nSym)

                IF (verbose) THEN
                        WRITE(out,'(a)') 'Atomic positions have been symmetrized'
                        CALL PrintSym3Dgroup(group, out)
                        WRITE(out,*)
                        WRITE(out,'(a)') '==========================='
                        WRITE(out,*)
                END IF
        END SUBROUTINE Symmetrize3DPositions

END MODULE Symmetry3DModule
