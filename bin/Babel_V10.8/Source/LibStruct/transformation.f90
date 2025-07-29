MODULE Transformation

  PUBLIC :: duplicate_cell, rotate_cell, translate_cell, clip_atoms

  INTEGER, parameter, private :: verbosity_max=3
  ! Zero for distances (in A)
  REAL(kind(0.d0)), parameter, private :: Distance_Zero=1.d-4
  REAL(kind(0.d0)), parameter, private :: Distance_Zero2=Distance_Zero*Distance_Zero

CONTAINS
  
  SUBROUTINE rotate_cell(rot, xp, im, at, out, verbose)

    USE MathStruct
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    REAL(kind(0.d0)), dimension(:,:), intent(inout) ::xp
    INTEGER, intent(in) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(inout) :: at
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    ! Identity matrix
    REAL(kind(0.d0)), dimension(1:3,1:3), parameter :: matId=ReShape( &
        (/ 1.d0, 0.d0, 0.d0, &
           0.d0, 1.d0, 0.d0, &
           0.d0, 0.d0, 1.d0 /), (/ 3,3 /) )

    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_rot
    REAL(kind(0.d0)) :: det

    ! Check that rot is a rotation matrix
    inv_rot = Transpose( rot )
    IF (mat3Norm2( MatMul( inv_rot, rot ) - matId ).GT.1.d-6) THEN
            WRITE(0,'(a)') 'Rot(1:3,1:3) is not a unitary matrix'
            STOP '< rotate_cell >'
    END IF
    ! Check that it is direct
    det=Mat3Det(rot(:,:))
    IF (Abs(det-1.d0).GT.1.d-6) THEN
            WRITE(0,'(a)') 'Determinant of Rot(1:3,1:3) is not equal to 1'
            WRITE(0,'(a,g18.6)') '  det = ', det
            IF (Abs(det+1.d0).LE.1.d-6) THEN
                    WRITE(0,'(a)') 'This corresponds to a mirror symmetry'
                    WRITE(0,'(a)') 'WARNING < rotate_cell >'
            ELSE
                    STOP '< Rotate >'
            END IF
    END IF
    
    ! Rotate atom cartesian coordinates
    IF (im.GT.0) xp(:,1:im) = MatMul( rot(:,:), xp(:,1:im) )

    ! Rotate unit cell vectors 
    at(:,:) = MatMul( rot, at )

    ! Print new values
    IF (verbose) THEN
            WRITE(out,'(a)') 'Unit cell has been rotated'
            IF (im.GT.0) WRITE(out,'(a,3g14.6)') '  gravity center located in: ', &
                    SUM(xp(1:3,1:im),2)/dble(im)
            IF (mat3Norm2(at).GT.1.d-6) THEN
                    WRITE(out,'(a)') '  new unit cell vectors:'
                    WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
                    WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
                    WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
            END IF
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
    END IF

  END SUBROUTINE rotate_cell

  SUBROUTINE Translate_cell(u, xp, im, out, verbose)

     IMPLICIT NONE

     ! Translation vector
     REAL(kind(0.d0)), dimension(1:3), intent(in) :: u
     ! Atom coordinates: xp(1:3,:) 
     REAL(kind(0.d0)), dimension(:,:), intent(inout) ::xp
     ! Number of atoms in simulation box 
     INTEGER, intent(in) :: im
     ! Output unit
     INTEGER, intent(in) :: out
     LOGICAL, intent(in) :: verbose

     INTEGER :: i, n

     ! Translate atom cartesian coordinates
     DO i=1, im
        xp(1:3,i) = xp(1:3,i) + u(1:3)
     END DO

     ! Print new values
     IF (verbose) THEN
             WRITE(out,'(a)') 'Unit cell has been translated'
             IF (im.GT.0) WRITE(out,'(a,3g14.6)') '  gravity center located in: ', &
                     SUM(xp(1:3,1:im),2)/dble(im)
             WRITE(out,*)
             WRITE(out,'(a)') '==========================='
             WRITE(out,*)
     END IF

  END SUBROUTINE Translate_cell

  SUBROUTINE duplicate_cell(lat, xp, iTyp, im, at, out, verbose)

    IMPLICIT NONE

    ! Number of times the unit cell is duplicated lat(1), lat(2)... in each direction
    INTEGER, dimension(3), intent(in) :: lat
    ! Atom coordinates: xp(1:3,:) 
    REAL(kind(0.d0)), dimension(:,:), intent(inout) ::xp
    ! Atom type
    INTEGER, dimension(:), intent(inout) :: iTyp
    ! Number of atoms in simulation box 
    INTEGER, intent(inout) :: im
    ! Lattice vector coordinates (A): at(1:3,1), at(2:3,2), ...
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(inout) :: at
    ! Output unit
    INTEGER, intent(in) :: out
    ! Verbose mode
    LOGICAL, intent(in) :: verbose

    INTEGER :: imm, new_im, i, j, k, l, n, iCell
    REAL(kind(0.d0)), dimension(1:3) :: center
    REAL(kind(0.d0)) :: at_norm2

    ! Maximal number of atoms
    imm = Min( Size(iTyp, 1), Size(xp, 2) )

    IF (verbose) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Before duplication of unit cell'
            WRITE(out,'(a,i0)') '  number of atoms: ', im
            IF (im.GT.0) WRITE(out,'(a,3g14.6)') '  gravity center located in: ', &
                    SUM(xp(1:3,1:im),2)/dble(im)
            WRITE(out,'(a)') '  unit cell vectors:'
            WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
            WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
            WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
    END IF

    IF ( (lat(1).EQ.0).OR.(lat(2).EQ.0).OR.(lat(3).EQ.0) ) THEN
            WRITE(0,'(a)') 'Numer of duplicated unit cell lat(1:3) cannot be null'
            STOP '< Duplicate_cell >'
    END IF

    at_norm2=SUM( at(1:3,1:3)**2 )
    IF (at_norm2.LE.Distance_Zero2) THEN
            WRITE(0,'(a)') 'Basis vectors at(1:3,1:3) undefined'
            WRITE(0,'(a)') ' => could not duplicate unit cell'
            WRITE(0,'(a,g20.12)') '  at_norm2 = ', at_norm2
            WRITE(0,'(a,g20.12)') '  distance_zero2 = ', distance_zero2
            STOP '< Duplicate_cell >'
    END IF

    ! Number of atoms after duplication
    new_im=im*Product(Abs(lat(1:3)))
    IF (new_im.GT.imm) THEN
            WRITE(0,'(a)') 'Problem when replicating unit cell'
            WRITE(0,'(a,i0)') '  max. number of atoms (imm): ', imm
            WRITE(0,'(a,i0)') '  needed number of atoms: ', new_im
            STOP '< Duplicate >'
    END IF

    ! Duplicate unit cell
    DO l=1, Abs(lat(3))
       DO k=1, Abs(lat(2))
          DO j=1, Abs(lat(1))
             iCell = im*( (l-1)*Abs(lat(2)*lat(1)) + (k-1)*Abs(lat(1)) +j-1 )
             IF (iCell.EQ.0) Cycle
             center(1:3) = dble(j-1)*at(1:3,1)*Dble(sign(1,lat(1))) &
                + dble(k-1)*at(1:3,2)*Dble(sign(1,lat(2))) &
                + dble(l-1)*at(1:3,3)*Dble(sign(1,lat(3)))
             DO i=1, im
                n = i + iCell
                xp(1:3,n ) = xp(1:3,i) + center(1:3)
                ityp(n) = ityp(i)
             END DO
          END DO
       END DO
    END DO

    im=new_im

    ! Lattice vector coordinates of the simulation box
    at(1:3,1) = Abs(dble(lat(1)))*at(1:3,1)
    at(1:3,2) = Abs(dble(lat(2)))*at(1:3,2)
    at(1:3,3) = Abs(dble(lat(3)))*at(1:3,3)

    IF (verbose) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Unit cell has been duplicated'
            WRITE(out,'(a,i0)') '  number of atoms: ', im
            IF (im.GT.0) WRITE(out,'(a,3g14.6)') '  gravity center located in: ', &
                    SUM(xp(1:3,1:im),2)/dble(im)
            WRITE(out,'(a)') '  unit cell vectors:'
            WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
            WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
            WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
    END IF

  END SUBROUTINE duplicate_cell

  SUBROUTINE Clip_Atoms(xp, im, at, out, verbose)
    ! Clip atoms in unit cell applying periodic boundary conditions
    
    USE MathStruct
    IMPLICIT NONE

    ! Atom coordinates: xp(1:3,:) 
    REAL(kind(0.d0)), dimension(:,:), intent(inout) ::xp
    ! Number of atoms in simulation box 
    INTEGER, intent(in) :: im
    ! Lattice vector coordinates (A): at(1:3,1), at(2:3,2), ...
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: at
    ! Output unit
    INTEGER, intent(in) :: out
    ! Verbose mode
    LOGICAL, intent(in) :: verbose

    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
    REAL(kind(0.d0)), dimension(1:3,1:im) :: xc
    REAL(kind(0.d0)) :: at_norm2

    at_norm2=SUM( at(1:3,1:3)**2 )
    IF (at_norm2.LE.Distance_Zero2) THEN
            WRITE(0,'(a)') 'You need to define lattice vectors at(1:3,i) to be&
                        & able to clip atoms'
            STOP '< Clip_Atoms >'
    END IF
    CALL Mat3Inv(at,inv_at)

    ! Apply periodic boundary conditions to atom coordinates
    xc(:,1:im) = MatMul( inv_at(:,:), xp(:,1:im) )
    xc(:,1:im) = Modulo( xc(:,1:im), 1.d0 )
    xp(:,1:im) = MatMul( at(:,:), xc(:,1:im) )

    IF (verbose) THEN
            WRITE(out,*)
            WRITE(out,'(a)') 'Periodic boundary conditions have been applied to atom coordinates'
            WRITE(out,*)
            WRITE(out,'(a)') '==========================='
            WRITE(out,*)
    END IF

  END SUBROUTINE Clip_Atoms

END MODULE Transformation
