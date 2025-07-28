MODULE elasticity_ani

  INTEGER, parameter, private :: verbosity_max=7

CONTAINS

  FUNCTION Isotropic_CVoigt(CVoigt) RESULT(test)
    ! test=.true. if the elastic constants CVoigt are isotropic

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:6,1:6), intent(in) :: CVoigt
    LOGICAL :: test

    INTEGER :: i, j
    REAL(kind(0.d0)) :: C11, C12, C44, Bulk_Modulus, elastic_zero

    test=.TRUE.
    Bulk_Modulus = Sum(CVoigt(1:3,1:3))/9.d0
    elastic_Zero = 1.d-10*Bulk_Modulus
    C11=CVoigt(1,1) ; C12=CVoigt(1,2) ; C44=CVoigt(4,4)

    IF (Abs(C11-C12-2.d0*C44).GT.elastic_zero) test=.FALSE.

    DO i=1, 3
       IF (Abs(CVoigt(i,i)-C11).GT.elastic_zero) test=.FALSE.
       DO j=i+1, 3
          IF (Abs(CVoigt(i,j)-C12).GT.elastic_zero) test=.FALSE.
          IF (Abs(CVoigt(j,i)-C12).GT.elastic_zero) test=.FALSE.
       END DO
       DO j=4, 6
          IF (Abs(CVoigt(i,j)).GT.elastic_Zero) test=.FALSE.
          IF (Abs(CVoigt(j,i)).GT.elastic_Zero) test=.FALSE.
       END DO
    END DO
    DO i=4, 6
       IF (Abs(CVoigt(i,i)-C44).GT.elastic_zero) test=.FALSE.
       DO j=i+1, 6
          IF (Abs(CVoigt(i,j)).GT.elastic_Zero) test=.FALSE.
          IF (Abs(CVoigt(j,i)).GT.elastic_Zero) test=.FALSE.
       END DO
    END DO

  END FUNCTION Isotropic_CVoigt

  !====================================================================

  FUNCTION Isotropic_Elastic_Constants(elastic_C) RESULT(test)
    ! test=.true. if the elastic constants are isotropic
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3), intent(in) :: elastic_C
    LOGICAL :: test
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt

    CVoigt=Pack_CVoigt(elastic_C)
    test=Isotropic_CVoigt( CVoigt )

  END FUNCTION Isotropic_Elastic_Constants

  !====================================================================

  FUNCTION Unpack_CVoigt(CVoigt) RESULT(elastic_C)
     ! Transfrom elastic constant from Voigt notation CVoigt(1:6,1:6) to full
     ! matrix elastic_C(i,j,k,l)
     IMPLICIT NONE
     REAL(kind(0.d0)), dimension(1:6,1:6), intent(in) :: CVoigt
     REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: elastic_C

    INTEGER :: i, j, k, l, m, n

    DO i=1,3 ; DO j=1,3
       IF (i.EQ.j) THEN 
               m=i
       ELSE
               m=9-i-j
       END IF
       DO k=1,3 ;  DO l=1,3
          IF (k.EQ.l) THEN 
                  n=l
          ELSE
                  n=9-k-l
          END IF
          elastic_C(i,j,k,l) = CVoigt(m,n)
       END DO ; END DO
    END DO ; END DO

  END FUNCTION Unpack_CVoigt

  !====================================================================

  FUNCTION CVoigt_random(Cmax, out) RESULT(CVoigt)
    ! Generate random elastic constant with maximal amplitude Cmax
    !  i.e. 6x6 definite positive tensor
    USE babel_data, ONLY : verbosity
    IMPLICIT NONE
    INTEGER, intent(in) :: out
    REAL(kind(0.d0)), intent(in) :: Cmax
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt

    INTEGER :: i, j
    REAL(kind(0.d0)) :: rnd, inv_norm, scalar
    REAL(kind(0.d0)), dimension(1:6,1:6) :: rot, inv_rot, ev

    IF (Cmax.LE.0.d0) THEN
            WRITE(0,*) "Cannot generate random elastic constants with a negative&
                & or null maximum value"
            WRITE(0,'(a,g12.3)') 'Cmax = ', Cmax
            STOP '< CVoigt_random >'
    END IF

    ! Generate 6 positive eigen-values
    ev(:,:) = 0.d0
    DO i=1, 6
       CALL Random_Number(rnd)
       ev(i,i) = Cmax*rnd
    END DO

    ! Generate random 6x6 rotation matrix
    !   - initialization
    rot(:,:) = 0.d0
    DO i=1, 6
       rot(i,i) = 1.d0
    END DO
    DO i=1, 6
       !  - random vector i
       CALL Random_Number(rnd)
       rot(:,i) = rnd*rot(:,i)
       DO j=i+1, 6
          CALL Random_Number(rnd)
          rot(:,i) = rot(:,i) + rnd*rot(:,j)
       END DO
       !  - normalize vector i
       inv_norm = 1.d0/Sqrt( Sum( rot(:,i)**2 ) )
       rot(:,i) = inv_norm*rot(:,i)
       !  - orthogonalize and normalize remaining vectors
       DO j=i+1, 6
          scalar = Sum( rot(:,i)*rot(:,j) )
          rot(:,j) = rot(:,j) - scalar*rot(:,i)
          inv_norm = 1.d0/Sqrt( Sum( rot(:,j)**2 ) )
          rot(:,j) = inv_norm*rot(:,j)
       END DO
    END DO

    ! Inverse rotation matrix
    inv_rot = Transpose( rot )

    ! Random definite positive matrix
    CVoigt = MatMul( inv_rot, MatMul( ev, rot) )
    
    IF (verbosity.GE.verbosity_max) THEN
            WRITE(out,*)
            WRITE(out,'(a,g14.6)') 'Generate random elastic constants with maximum value ', Cmax
            WRITE(out,'(a,6g14.6)') '  eigen-values: ', (ev(i,i), i=1, 6)
            WRITE(out,'(a)') '  rotation matrix:'
            DO i=1, 6
               WRITE(out,'(a,6g14.6,a)') '     | ', rot(i,1:6), ' |'
            END DO
            WRITE(out,'(a)') '  random definite positive matrix:'
            DO i=1, 6
               WRITE(out,'(a,6g14.6,a)') '     | ', CVoigt(i,1:6), ' |'
            END DO
    END IF


  END FUNCTION CVoigt_random

  !====================================================================

  SUBROUTINE Check_CVoigt(CVoigt) 
    ! Check if matrix CVoigt is symmetric and definite positive

    USE Math
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:6,1:6), intent(in) :: CVoigt

    REAL(kind(0.d0)), dimension(1:6,1:6) :: A
    REAL(kind(0.d0)), dimension(1:6) :: W
    REAL(kind(0.d0)), dimension(1:17) :: Work
    INTEGER :: info

    ! Symmetric matrix
    IF (MatNorm2( CVoigt - Transpose(CVoigt)).GE.1.d-6) THEN
            WRITE(0,'(a)') 'Elastic constants CVoigt(1:6,1:6) need to &
            &be symmetric'
            CALL Print_CVoigt(CVoigt,0)
            STOP '< Check_CVoigt >'
    END IF

    ! Calculate matrix eigenvalues using dsyev subroutine from Lapack package
    A(:,:) = CVoigt(:,:)
    CALL Dsyev('N', 'U', 6, A, 6, W, work, 17, info)
    IF (Info.NE.0) THEN
            WRITE(0,'(a)') 'Error when calling Dsyev subroutine from Lapack package'
            WRITE(0,'(a,i0)') '  info = ', info
    END IF

    ! Check that smaller eigenvalue is positive
    IF (W(1).LE.0.d0) THEN
            WRITE(0,'(a)') 'Matrix defining elastic constants is not definite positive'
            CALL Print_CVoigt(CVoigt,0)
            WRITE(0,'(a,6g14.6)') ' eigenvalues: ', W(1:6)
            STOP '< Check_CVoigt >'
    END IF

  END SUBROUTINE Check_CVoigt

  !====================================================================

  FUNCTION Pack_CVoigt(elastic_C) RESULT(CVoigt)
    ! Transfrom elastic constant from full
    ! matrix elastic_C(i,j,k,l) to Voigt notation CVoigt(1:6,1:6) 
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3), intent(in) :: elastic_C
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt

    INTEGER :: i, j, k, l, m, n

    DO i=1,3 ; DO j=i,3
       IF (i.EQ.j) THEN 
               m=i
       ELSE
               m=9-i-j
       END IF
       DO k=1,3 ;  DO l=k,3
          IF (k.EQ.l) THEN 
                  n=l
          ELSE
                  n=9-k-l
          END IF
          CVoigt(m,n) = 0.25d0*( elastic_C(i,j,k,l) + elastic_C(j,i,k,l) &
                + elastic_C(i,j,l,k) + elastic_C(j,i,l,k) )
       END DO ; END DO
    END DO ; END DO

  END FUNCTION Pack_CVoigt

  !====================================================================

  SUBROUTINE Print_CVoigt(CVoigt, out)
    ! Print the  6x6 Voigt matrix of elastic constants
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:6,1:6), intent(in) :: CVoigt
    INTEGER, intent(in) :: out  ! output unit

    INTEGER :: m

    DO m=1, 6
       WRITE(out,'(a,6(1x,g12.5),a)') '   | ', CVoigt(m,1:6), ' |'
    END DO

  END SUBROUTINE Print_CVoigt

  !====================================================================

  SUBROUTINE Print_Elastic_Constants(elastic_c, out)
    ! Print the  6x6 Voigt matrix of elastic constants
    IMPLICIT NONE
    INTEGER, intent(in) :: out  ! output unit
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3), intent(in) :: elastic_C
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt

    CVoigt=Pack_CVoigt(elastic_C)
    CALL Print_CVoigt( CVoigt, out )

  END SUBROUTINE Print_Elastic_Constants

  !====================================================================

  FUNCTION Rotate_CVoigt(CVoigt1, rot) RESULT(CVoigt2)
    ! Rotate elastic constants CVoigt(1:6,1:6) 
    IMPLICIT NONE
    ! Inverse of the rotation matrix
    REAL(kind(0.d0)), dimension(1:6,1:6), intent(in) :: CVoigt1
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    REAL(kind(0.d0)), dimension(1:6,1:6) :: CVoigt2

    ! Elastic constants in cartesian coordinates
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: elastic_C1, elastic_C2 

    elastic_C1 = Unpack_CVoigt( CVoigt1 )
    elastic_C2 = rotate_elastic_constants(elastic_C1, rot)
    CVoigt2 = Pack_CVoigt( elastic_C2 )

  END FUNCTION Rotate_CVoigt

  !====================================================================

  FUNCTION Rotate_Elastic_Constants(elastic_C1, rot) RESULT(elastic_C2)
    ! Rotate elastic constants elastic_C(1:3,1:3,1:3,1:3)
    IMPLICIT NONE
    ! Inverse of the rotation matrix
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3), intent(in) :: elastic_C1
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: elastic_C2 

    ! Elastic constants in cartesian coordinates
    INTEGER :: i, j, k, l, m, n, o, p

    ! Elastic matrix in the dislo axes
    !   Hirth and Lothe, pp.431-432
    elastic_C2=0.d0
    DO l=1, 3 ; DO k=1, 3 ; DO j=1, 3 ; DO i=1, 3
           DO p=1, 3 ; DO o=1, 3 ; DO n=1, 3 ; DO m=1, 3
                  elastic_C2(i,j,k,l) = elastic_C2(i,j,k,l) + &
                        rot(i,m)*rot(j,n)*rot(k,o)*rot(l,p)*elastic_C1(m,n,o,p)
           END DO ; END DO ; END DO ; END DO
    END DO ; END DO ; END DO ; END DO

  END FUNCTION Rotate_Elastic_Constants

  !====================================================================

  FUNCTION Rotate_eVoigt(eVoigt1, rot) RESULT(eVoigt2)
    ! Rotate strain keeping Voigt notation
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:6), intent(in) :: eVoigt1
    REAL(kind(0.d0)), dimension(1:3,1:3), intent(in) :: rot
    REAL(kind(0.d0)), dimension(1:6) :: eVoigt2

    REAL(kind(0.d0)), dimension(1:3,1:3) :: epsi1, epsi2, inv_rot

    inv_rot=Transpose(rot)
    epsi1(:,:) = Reshape( (/ eVoigt1(1), 0.5d0*eVoigt1(6), 0.5d0*eVoigt1(5), &
                             0.5d0*eVoigt1(6), eVoigt1(2), 0.5d0*eVoigt1(4), &
                             0.5d0*eVoigt1(5), 0.5d0*eVoigt1(4), eVoigt1(3) /), (/3,3/) )
    epsi2(:,:) = MatMul( rot(:,:), MatMul( epsi1(:,:), inv_rot(:,:) ) )
    eVoigt2(1) = epsi2(1,1)
    eVoigt2(2) = epsi2(2,2)
    eVoigt2(3) = epsi2(3,3)
    eVoigt2(4) = epsi2(2,3)+epsi2(3,2)
    eVoigt2(5) = epsi2(1,3)+epsi2(3,1)
    eVoigt2(6) = epsi2(1,2)+epsi2(2,1)

  END FUNCTION Rotate_eVoigt

  !====================================================================

  SUBROUTINE inv_Voigt_index(m,i,j)
    ! input: m, Voigt index in [1:9]
    ! outputs: i and j, corresponding 2D indexes in [1:3]
    IMPLICIT NONE
    INTEGER, intent(in) :: m
    INTEGER, intent(out) :: i, j

    SELECT CASE(m)
    CASE(1:3)
            i=m ; j=m
    CASE(4)
            i=2 ; j=3
    CASE(5)
            i=1 ; j=3
    CASE(6)
            i=1 ; j=2
    CASE(7)
            i=3 ; j=2
    CASE(8)
            i=3 ; j=1
    CASE(9)
            i=2 ; j=1
    END SELECT
  END SUBROUTINE inv_Voigt_index

END MODULE elasticity_ani
