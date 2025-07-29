MODULE SlabModule

  IMPLICIT NONE
  LOGICAL, save :: slab
  INTEGER, save :: nxSlab, nySlab, nzSlab

  REAL(kind(0.d0)), parameter, private :: Distance_Zero=1.d-4


CONTAINS

  SUBROUTINE MakeSlab(xp, iTyp, im, at, out, verbose)

    USE Math
    IMPLICIT NONE

    REAL(kind(0.d0)), dimension(:,:), intent(inout) ::xp
    INTEGER, dimension(:), intent(inout) :: iTyp
    INTEGER, intent(inout) :: im
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(inout) :: at
    INTEGER, intent(in) :: out
    LOGICAL, intent(in) :: verbose

    REAL(kind(0.d0)), dimension(:,:), allocatable :: xc
    REAL(kind(0.d0)), dimension(1:3,1:3) :: inv_at
    REAL(kind(0.d0)), dimension(1:3) :: uShift, sShift
    REAL(kind(0.d0)) :: inv_nxSlab, inv_nySlab, inv_nzSlab
    INTEGER :: i, im_new

    ! Check if periodicity vectors are defined
    IF (matNorm2(at(1:3,1:3)).LE.Distance_Zero) THEN
            WRITE(0,'(a)') 'Periodicity vectors are not defined'
            STOP '< MakeSlab >'
    END IF

    ! Calculate atom reduced coordinates
    ALLOCATE(xc(1:3,1:im))
    CALL MatInv(at,inv_at)
    xc(:,1:im) = MatMul( inv_at(:,:), xp(:,1:im) )
    xc(:,1:im) = Modulo( xc(:,1:im), 1.d0 )

    ! Add a small shift to all atoms
    uShift(:) = distance_zero * (/ 1.d0, 1.d0, 1.d0 /)
    sShift(:) = MatMul( inv_at(:,:), uShift(:) )
    DO i=1, im
       xc(:,i) = xc(:,i) + sShift(:)
    END DO

    ! Keep only atoms for which 0 < xc(1,:) < 1/nxSlab, ...
    inv_nxSlab = 1.d0/dble(nxSlab)
    inv_nySlab = 1.d0/dble(nySlab)
    inv_nzSlab = 1.d0/dble(nzSlab)

    im_new=0
    DO i=1, im
       IF ( (xc(1,i).LT.inv_nxSlab).AND.(xc(2,i).LT.inv_nySlab) &
                .AND.(xc(3,i).LT.inv_nzSlab) ) THEN
                ! keep atom
               im_new = im_new +1
               xp(:,im_new) = xp(:,i)
               iTyp(im_new) = iTyp(i)
       END IF
    END DO
    DEALLOCATE(xc)
    xp(:,im_new+1:im) = 0.d0
    iTyp(im_new+1:im) = 0

    im = im_new

    ! New periodicity vectors
    at(:,1) = inv_nxSlab*at(:,1)
    at(:,2) = inv_nySlab*at(:,2)
    at(:,3) = inv_nzSlab*at(:,3)

    IF (.NOT.verbose) RETURN
    WRITE(out,*)
    WRITE(out,'(a)') 'A slab has been extracted from unit cell'
    WRITE(out,'(a,3(i0,a))') '  periodicity vectors have been divided by ' , &
        nxSlab, ', ', nySlab, ', and ', nzSlab, ' respectively'
    WRITE(out,'(a,i0)') '  new number of atoms: ', im
    IF (im.GT.0) WRITE(out,'(a,3g14.6)') '  gravity center now located in: ', &
            SUM(xp(1:3,1:im),2)/dble(im)
    WRITE(out,'(a)') '  new unit cell vectors:'
    WRITE(out,'(a,3g20.12)') '    at(1:3,1): ', at(1:3,1)
    WRITE(out,'(a,3g20.12)') '    at(1:3,2): ', at(1:3,2)
    WRITE(out,'(a,3g20.12)') '    at(1:3,3): ', at(1:3,3)
    WRITE(out,*)
    WRITE(out,'(a)') '==========================='
    WRITE(out,*)

  END SUBROUTINE MakeSlab

END MODULE SlabModule
