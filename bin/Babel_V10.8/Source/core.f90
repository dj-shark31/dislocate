MODULE core_module
  ! This is not working with periodic boundary conditions

  INTEGER, parameter, private :: verbosity_max=4

CONTAINS

  SUBROUTINE Core(coreAtom, out)
    ! Determines atom that are in the core of line defects
    USE Babel_data
    USE disloc_elasticity_ani 
    USE lineCouple_elasticity_ani
    USE lineForce_elasticity_ani
    IMPLICIT NONE
    LOGICAL, dimension(:), intent(out) :: coreAtom
    INTEGER, intent(in) :: out

    INTEGER :: n, i
    REAL(kind(0.d0)) :: rc2
    REAL(kind(0.d0)), dimension(1:3) :: dR

    coreAtom(:)=.FALSE.

    IF (rc.LE.0.d0) RETURN
    rc2=rc**2

    DO n=1, nd
      DO i=1, im
        dR(:) = xp(:,i) - cDislo(:,n)
        dR(:) = dR(:) - Sum( dR(:)*lDislo(:,n) )*lDislo(:,n)
        IF ( Sum(dR(:)**2).LE.rc2 ) coreAtom(i)=.TRUE.
      END DO
    END DO

    DO n=1, nlf
      DO i=1, im
        dR(:) = xp(:,i) - cLineForce(:,n)
        dR(:) = dR(:) - Sum( dR(:)*lLineForce(:,n) )*lLineForce(:,n)
        IF ( Sum(dR(:)**2).LE.rc2 ) coreAtom(i)=.TRUE.
      END DO
    END DO

    DO n=1, nlc
      DO i=1, im
        dR(:) = xp(:,i) - cLineCouple(:,n)
        dR(:) = dR(:) - Sum( dR(:)*lLineCouple(:,n) )*lLineCouple(:,n)
        IF ( Sum(dR(:)**2).LE.rc2 ) coreAtom(i)=.TRUE.
      END DO
    END DO

    IF (verbosity.LT.verbosity_max) RETURN
    WRITE(out,*)
    WRITE(out,'(a)') '==========================='
    WRITE(out,*)
    WRITE(out,'(a,g14.6)') 'Radius cutoff used to define the core of line defects: rc = ', rc
    WRITE(out,'(a,i0)') '  number of atoms in cores: ', Count(coreAtom(:))
    WRITE(out,*)

  END SUBROUTINE Core

END MODULE core_module
