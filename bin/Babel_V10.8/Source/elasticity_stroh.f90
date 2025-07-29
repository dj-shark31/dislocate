MODULE elasticity_Stroh
  !------------------------------------------------------------------------------------
  ! Anisotropic elasticity for the case where the elastic state is independent
  ! ot the z Cartesian coordinates
  !     => used for dislocations, line-forces and line-force couples
  !----------------------------------------
  !   Ref.: J.P. Hirth and J. Lothe, Theory of Dislocations 
  !             (Wiley, New York, 2nd ed., 1982), pp.423-483
  !         J.P. Hirth and J. Lothe, 
  !             Anisotropic elastic solutions for line defects in high-symmetry cases
  !             J. Appl. Phys. 44 (1973), p. 1029
  !         J.D. Eshelby, W.T. Read and W. Shockley,
  !             Anisotropic elasticity with applications to dislocation theory
  !             Acta Met. 1 (1953), p. 251
  !         A.N. Stroh
  !             Dislocations and cracks in anisotropic elasticity
  !             Phil. Mag. 3 (1958), p. 625
  !         A.N. Stroh
  !             Steady state problems in anisotropic elasticity
  !             J. Math. Phys. 41 (1962), p. 77
  !   USE lapack subroutines: zgeev

  PUBLIC :: InitStroh, Build_DStroh
  
  !---------------------------------------------------------
  REAL(kind(0.d0)), save, private :: elastic_Zero

  INTEGER, parameter, private :: verbosity_max=6

CONTAINS

  !====================================================================

  SUBROUTINE InitStroh(elastic_Cp, root, elastic_A, elastic_B, elastic_L, out)
     ! - build the a(i,j) matrix of Eq. 13.57 p. 438 
     ! - solve the equation det[a(i,j)]=0
     ! - determine vector solution of a(i,j)*v(j)=0. => elastic_A
     ! - build the b(i,j,k) matrix of Eq. 13.77 p. 443 => elastic_B

     USE babel_data, ONLY : verbosity
     USE elasticity_ani
     USE math
     IMPLICIT NONE

     ! Elastic constants in dislo or line-force axes
     REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3), intent(in) :: elastic_Cp

     ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
     !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
     !   and the equation det[a(i;j)]=0 has 6 solutions for the variable p: root(1:6) 
     COMPLEX(kind(0.d0)), dimension(1:6), intent(out) :: root
     !   corresponding vector solution of the equation a(i,j)*elastic_A(j) = 0.
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(out) :: elastic_A
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(out) :: elastic_L
     ! Matrix B(i,j,k) of Eq. 13.77 p. 443 given by
     !   b(i,j,k) = Cijk1 + Cijk2*p
     COMPLEX(kind(0.d0)), dimension(1:3,1:3,1:3,1:6), intent(out) :: elastic_B

     INTEGER, intent(in), optional :: out
     INTEGER :: i,j, m,n,p
     REAL(kind(0.d0)) :: factor, det_norm2
     COMPLEX(kind(0.d0)) :: cfactor
     ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
     !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
     COMPLEX(kind(0.d0)), dimension(1:3,1:3,1:6) :: a
     REAL(kind(0.d0)), dimension(0:2,1:3,1:3) :: a_poly
     ! Coef of the polynomial corresponding to the determinant of the a(i,j) matrix
     REAL(kind(0.d0)), dimension(0:6) :: det_poly
     ! Companion matrix of this polynomial
     COMPLEX(kind(0.d0)), dimension(1:6,1:6) :: det_poly_companion
     ! Variables used by lapack subroutine
     CHARACTER(len=1) :: Jobvl, Jobvr
     INTEGER :: lda, ldvl, ldvr, lwork, info
     COMPLEX(kind(0.d0)), dimension(1:12) :: WORK
     REAL(kind(0.d0)), dimension(1:12) :: RWORK
     COMPLEX(kind(0.d0)), dimension(6,1) :: VL, VR
     COMPLEX(kind(0.d0)), dimension(6) :: eigenValue
     REAL(kind(0.d0)), dimension(1:6) :: RWORK2
     COMPLEX(kind(0.d0)), dimension(1:3,1:3) :: atemp
     COMPLEX(kind(0.d0)), dimension(1:3) :: aev
     COMPLEX(kind(0.d0)), dimension(1:3) :: vector
     COMPLEX(kind(0.d0)) :: det, inv_norm
     COMPLEX(kind(0.d0)), dimension(1:6,1:6) :: matAL
     COMPLEX(kind(0.d0)), dimension(1:3,1:3) :: vecAL, vecAA, vecLL
     LOGICAL :: test
     REAL(kind(0.d0)) :: Bulk_Modulus

     IF (isotropic_elastic_constants(elastic_Cp)) THEN
             WRITE(0,'(a)') 'Sextic formalism cannot handle isotropic elasticity'
             STOP '< InitStroh >'
     END IF

     Bulk_Modulus=0.d0
     DO i=1, 3 ;  DO j=1, 3
        Bulk_Modulus = Bulk_Modulus + elastic_Cp(i,i, j,j)
     END DO ;  END DO
     Bulk_Modulus = Bulk_Modulus/9.d0
     elastic_Zero = 1.d-12*Bulk_Modulus

     IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
             WRITE(out,*)
             WRITE(out,'(a)') 'Initialization for line defect using Stroh formalism'
             WRITE(out,'(a)') '  elastic constants in line axes'
             CALL Print_Elastic_Constants(elastic_Cp,out)
     END IF

     ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
     !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
     DO i=1,3 ; DO j=1,3
        a_poly(0,i,j) = elastic_Cp(i,1,j,1)
        a_poly(1,i,j) = elastic_Cp(i,1,j,2) + elastic_Cp(i,2,j,1)
        a_poly(2,i,j) = elastic_Cp(i,2,j,2)
     END DO ; END DO

     
     ! Coef of the polynomial corresponding to the determinant of the a(i,j) matrix
     det_poly(0:6) = 0.d0
     DO m=0,2 ; DO n=0,2 ; DO p=0,2
        det_poly(m+n+p) = det_poly(m+n+p) &
                + a_poly(m,1,1)*a_poly(n,2,2)*a_poly(p,3,3) &
                - a_poly(m,1,1)*a_poly(n,2,3)*a_poly(p,3,2) &
                - a_poly(m,1,2)*a_poly(n,2,1)*a_poly(p,3,3) &
                + a_poly(m,1,2)*a_poly(n,3,1)*a_poly(p,2,3) &
                + a_poly(m,1,3)*a_poly(n,2,1)*a_poly(p,3,2) &
                - a_poly(m,1,3)*a_poly(n,3,1)*a_poly(p,2,2) 
     END DO ; END DO ; END DO

     ! Monic polonomial
     factor = 1.d0/det_poly(6)
     det_poly(0:6) = factor*det_poly(0:6)

     ! Companion matrix of the polynomial
     !   => its eigenvalues are the roots of the polynomial
     det_poly_companion(1,1:6) = (/ (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), &
        (0.d0,0.d0), (0.d0,0.d0), Cmplx(-det_poly(0),0.d0,kind(0.d0)) /)
     det_poly_companion(2,1:6) = (/ (1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), &
        (0.d0,0.d0), (0.d0,0.d0), Cmplx(-det_poly(1),0.d0,kind(0.d0)) /)
     det_poly_companion(3,1:6) = (/ (0.d0,0.d0), (1.d0,0.d0), (0.d0,0.d0), &
        (0.d0,0.d0), (0.d0,0.d0), Cmplx(-det_poly(2),0.d0,kind(0.d0)) /)
     det_poly_companion(4,1:6) = (/ (0.d0,0.d0), (0.d0,0.d0), (1.d0,0.d0), &
        (0.d0,0.d0), (0.d0,0.d0), Cmplx(-det_poly(3),0.d0,kind(0.d0)) /)
     det_poly_companion(5,1:6) = (/ (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), &
        (1.d0,0.d0), (0.d0,0.d0), Cmplx(-det_poly(4),0.d0,kind(0.d0)) /)
     det_poly_companion(6,1:6) = (/ (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), &
        (0.d0,0.d0), (1.d0,0.d0), Cmplx(-det_poly(5),0.d0,kind(0.d0)) /)

     ! Call Lapack subroutine to find eigenvalues of the matrix det_poly_companion
     Jobvl='N'  ! Left eigenvectors are not computed
     Jobvr='N'  ! Right eigenvectors are not computed
     n=6        ! Order of the matrix det_poly_companion
     lda=6      ! Leading dimension of the array det_poly_companion
     ! eigenValue contains the computed eigenvalues
     ! Vl and Vr not used
     ldvl=1     ! leading dimension of the array vl
     ldvr=1     ! leading dimension of the array vr
     lwork=12   ! leading dimension of the array work (>=max(1,2*n))
     CALL Zgeev(JOBVL, JOBVR, N, det_poly_companion, LDA, eigenValue, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )

     ! Sort roots according to the rules: 
     !     - root(1)=root(4)*, root(2)=root(5)*, root(3)=root(6)*
     !     where * means complex conjugate
     !     - root(1), root(2), and root(3) have positive imaginary part
     root(1) = Cmplx(dble(eigenValue(1)),abs(aimag(eigenValue(1))), kind(0.d0))
     root(4) = Cmplx(dble(eigenValue(1)),-abs(aimag(eigenValue(1))), kind(0.d0))
     IF ( ( iNormInf( eigenValue(2) - root(1) ).LE.elastic_Zero ).OR. &
        ( iNormInf( eigenValue(2) - root(4) ).LE.elastic_Zero ) ) THEN
             ! eigenValue(2) is equal to root(1) or root(4)
             n=3
     ELSE
             n=2
     END IF
     root(2) = Cmplx(dble(eigenValue(n)),abs(aimag(eigenValue(n))), kind(0.d0))
     root(5) = Cmplx(dble(eigenValue(n)),-abs(aimag(eigenValue(n))), kind(0.d0))
     n=n+1
     test=.TRUE.
     DO WHILE (test)
         IF ( ( iNormInf( eigenValue(n) - root(1) ).LE.elastic_Zero ).OR. &
            ( iNormInf( eigenValue(n) - root(4) ).LE.elastic_Zero ).OR. &
            ( iNormInf( eigenValue(n) - root(2) ).LE.elastic_Zero ).OR. &
            ( iNormInf( eigenValue(n) - root(5) ).LE.elastic_Zero ) ) THEN
                 ! eigenValue(n) is equal to root(1) or root(4) or root(2) or root(5)
                 n=n+1
         ELSE
                 test=.FALSE.
         END IF
     END DO
     root(3) = Cmplx(dble(eigenValue(n)),abs(aimag(eigenValue(n))), kind(0.d0))
     root(6) = Cmplx(dble(eigenValue(n)),-abs(aimag(eigenValue(n))), kind(0.d0))

     DO n=1, 6
        ! Build the matrix a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
        DO j=1, 3 ; DO i=1, 3
           a(i,j,n) = a_poly(0,i,j) + ( a_poly(1,i,j) + a_poly(2,i,j)*root(n) )*root(n)
        ENDDO ; ENDDO
        ! Vectors of the matrix 
        ! for the six different values of the variable p(1:6)=roots(1:6)
        elastic_A(1,n) = a(1,2,n)*a(2,3,n) - a(2,2,n)*a(1,3,n)
        elastic_A(2,n) = - a(1,1,n)*a(2,3,n) + a(2,1,n)*a(1,3,n)
        elastic_A(3,n) = a(1,1,n)*a(2,2,n) - a(2,1,n)*a(1,2,n)
     END DO

     DO n=1, 6  ! Loop on the six roots 
        ! Matrix B(i,j,k) of Eq. 13.77 p. 443 given by
        !   b(i,j,k) = Cijk1 + Cijk2*p
        elastic_B(1:3,1:3,1:3, n) = elastic_Cp(1:3,1:3,1:3,1) + elastic_Cp(1:3,1:3,1:3,2)*root(n)

        ! Matrix L defined by Stroh
        elastic_L(1:3,n) = MatMul( elastic_B(1:3,2,1:3,n), elastic_A(1:3,n) )

        ! Normalize elastic_A and elastic_L according to Stroh orthogonality relation 2A.L = 1
        inv_norm = 1.d0/Sqrt( 2.d0*Sum( elastic_A(1:3,n) * elastic_L(1:3,n) ) )
        elastic_A(1:3,n) = inv_norm*elastic_A(1:3,n)
        elastic_L(1:3,n) = inv_norm*elastic_L(1:3,n)

     ENDDO      ! Loop on the six roots 

     IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
             WRITE(out,*)
             WRITE(out,'(a)') 'Roots before sorting: '
             DO n=1, 6
                WRITE(out,'(a,i1,2(a,g14.6))') '  n=', n, ' :  ', &
                        Dble(eigenvalue(n)), ' + i*', Aimag(eigenvalue(n))
             END DO
             WRITE(out,'(a)') 'Roots after sorting: '
             DO n=1, 6
                WRITE(out,'(a,i1,2(a,g14.6))') '  n=', n, ' :  ', &
                     Dble(root(n)), ' + i*', Aimag(root(n))
             END DO
             DO n=1, 6
                WRITE(out,*)
                WRITE(out,'(a,i1,a,g13.6,a,g13.6)') 'Root ', n, ': ', Dble(root(n)), ' + i*', Aimag(root(n))                
                WRITE(out,'(a)') '  corresponding a matrix'
                DO i=1,3
                   WRITE(out,'(a,3(g11.4,a,g11.4,2x),a)') &
                        '    | ', Dble(a(i,1,n)) , ' + i*', Aimag(a(i,1,n)), &
                        Dble(a(i,2,n)) , ' + i*', Aimag(a(i,2,n)), &
                        Dble(a(i,3,n)) , ' + i*', Aimag(a(i,3,n)), ' |'
                END DO
                det=matdet(a(:,:,n))
                WRITE(out,'(2(a,g14.6))') '  determinant: ', Dble(det), ' + i*', Aimag(det)
                ! Calculate eigen values of A matrix
                atemp(:,:) = a(:,:,n)
                CALL Zgeev('N', 'N', 3, atemp, 3, aev, VL, 1, VR, 1, WORK, 12, RWORK2, INFO )
                WRITE(out,'(a,6(g11.4,a))') '  corresponding eigenvalues : ( ', &
                        Dble(aev(1)), ' + i*', Aimag(aev(1)), ' , ', &
                        Dble(aev(2)), ' + i*', Aimag(aev(2)), ' , ', &
                        Dble(aev(3)), ' + i*', Aimag(aev(3)), ' )'
                WRITE(out,'(a,6(g11.4,a))') '  eigenvector of null eigenvalue: A = ( ', &
                        Dble(elastic_A(1,n)), ' + i*', Aimag(elastic_A(1,n)), ' , ', &
                        Dble(elastic_A(2,n)), ' + i*', Aimag(elastic_A(2,n)), ' , ', &
                        Dble(elastic_A(3,n)), ' + i*', Aimag(elastic_A(3,n)), ' )'
                ! Check eigenVector
                Vector(1:3) = MatMul(a(1:3,1:3,n), elastic_A(1:3,n)) 
                WRITE(out,'(a,6(g11.4,a))') '  a * A =  ( ', &
                        Dble(Vector(1)), ' + i*', Aimag(Vector(1)), ' , ', &
                        Dble(Vector(2)), ' + i*', Aimag(Vector(2)), ' , ', &
                        Dble(Vector(3)), ' + i*', Aimag(Vector(3)), ' )'
                WRITE(out,'(a,6(g11.4,a))') '  associated vector L = ( ', &
                        Dble(elastic_L(1,n)), ' + i*', Aimag(elastic_L(1,n)), ' , ', &
                        Dble(elastic_L(2,n)), ' + i*', Aimag(elastic_L(2,n)), ' , ', &
                        Dble(elastic_L(3,n)), ' + i*', Aimag(elastic_L(3,n)), ' )'
             END DO
             ! Check orthogonality relation between A and L
             DO n=1, 6 ; DO m=1, 6
                matAl(n,m) = Sum( elastic_A(1:3,n)*elastic_L(1:3,m) )
             END DO ; END DO
             matAL(1:6,1:6) = matAL(1:6,1:6) + Transpose(matAL(1:6,1:6))
             WRITE(out,*)
             WRITE(out,'(a)') 'Orthogonaly between A and L vectors'
             WRITE(out,*)
             WRITE(out,'(a)') "  Sum_i { A(i,a)*L(i,b) + A(i,b)*L(i,a) } :"
             DO n=1,6
                   WRITE(out,'(a,6(g9.2,a,g9.2,2x),a)') &
                        '    | ', (Dble(matAL(n,m)) , ' + i*', Aimag(matAL(n,m)), m=1,6), ' |'
             END DO

             WRITE(out,*)
             WRITE(out,'(a)') "  Sum_a { A(i,a)*A(j,a) } :"
             DO i=1, 3
                DO j=1, 3
                   vecAA(i,j) = Sum( elastic_A(i,1:6)*elastic_A(j,1:6) )
                END DO
                WRITE(out,'(a,3(g9.2,a,g9.2,2x),a)') &
                        '    | ', (Dble(vecAA(i,j)) , ' + i*', Aimag(vecAA(i,j)), j=1,3), ' |'
             END DO
             WRITE(out,*)
             WRITE(out,'(a)') "  Sum_a { L(i,a)*L(j,a) } :"
             DO i=1, 3
                DO j=1, 3
                   vecLL(i,j) = Sum( elastic_L(i,1:6)*elastic_L(j,1:6) )
                END DO
                WRITE(out,'(a,3(g9.2,a,g9.2,2x),a)') &
                        '    | ', (Dble(vecLL(i,j)) , ' + i*', Aimag(vecLL(i,j)), j=1,3), ' |'
             END DO
             WRITE(out,*)
             WRITE(out,'(a)') "  Sum_a { A(i,a)*L(j,a) } :"
             DO i=1, 3
                DO j=1, 3
                   vecAL(i,j) = Sum( elastic_A(i,1:6)*elastic_L(j,1:6) )
                END DO
                WRITE(out,'(a,3(g9.2,a,g9.2,2x),a)') &
                        '    | ', (Dble(vecAL(i,j)) , ' + i*', Aimag(vecAL(i,j)), j=1,3), ' |'
             END DO

             WRITE(out,*)
             WRITE(out,'(a)') 'Matrix appearing in energy prefactor for couples of line-force'
             DO n=1, 6 ; DO m=1, 6
                matAl(n,m) = Sum( elastic_A(1:3,n)*elastic_L(1:3,m) )
             END DO ; END DO
             matAL(1:6,1:6) = matAL(1:6,1:6) - Transpose(matAL(1:6,1:6))
             WRITE(out,'(a)') "  Sum_i { A(i,a)*L(i,b) - A(i,b)*L(i,a) } :"
             DO n=1,6
                   WRITE(out,'(a,6(g9.2,a,g9.2,2x),a)') &
                        '    | ', (Dble(matAL(n,m)) , ' + i*', Aimag(matAL(n,m)), m=1,6), ' |'
             END DO
             WRITE(out,*)

             WRITE(out,'(a)') 'Matrix appearing in dislocation interaction'
             DO n=1, 6
                WRITE(out,'(a,i0)') ' L(i,a)*L(j,a), a = ', n
                DO i=1, 3
                   WRITE(out,'(a,3(g9.2,a,g9.2,2x),a)') &
                        '    | ', (Dble(elastic_L(i,n)*elastic_L(j,n)) , & 
                        ' + i*', Aimag(elastic_L(i,n)*elastic_L(j,n)), j=1,3), ' |'
                END DO
             END DO

     END IF

     ! Check that vectors A are linearly independant
     det=MatDet( elastic_A(1:3,1:3) )
     det_norm2 = Sqrt( Dble(det)**2 + aImag(det)**2 )
     IF ( det_norm2 .LE. 1.d-7*matNorm2(elastic_A(1:3,1:3)) ) THEN
             WRITE(0,*)
             WRITE(0,'(a)') 'WARNING < InitStroh >'
             WRITE(0,'(a)') 'Vectors A(i,a) are not linearly independant'
             WRITE(0,'(a,g14.6)') '| Det(A) |^2 = ', det_norm2
             WRITE(0,'(a,g14.6)') '| A |^2 = ', matNorm2(elastic_A(1:3,1:3))
             WRITE(0,'(a)') 'You need to add some noise to elastic &
                &constants so as to break symmetry'
             WRITE(0,'(a)') 'Add to your input file "CVoigt_noise=1.d-4"'
            WRITE(0,'(a)') '  (do not use a noise lower than 1.d-6)'
             WRITE(0,*)
             STOP '< InitStroh >'
     ELSE IF (verbosity.GE.verbosity_max) THEN
             WRITE(out,'(a)') 'Vectors A(i,a) are linearly independant'
     END IF

     ! Check that vectors L are linearly independant
     det=MatDet( elastic_L(1:3,1:3) )
     det_norm2 = Sqrt( Dble(det)**2 + aImag(det)**2 )
     IF ( det_norm2 .LE. 1.d-7*matNorm2(elastic_L(1:3,1:3)) ) THEN
             WRITE(0,*)
             WRITE(0,'(a)') 'WARNING < InitStroh >'
             WRITE(0,'(a)') 'Vectors L(i,a) are not linearly independant'
             WRITE(0,'(a,g14.6)') '| Det(L) |^2 = ', det_norm2
             WRITE(0,'(a,g14.6)') '| L |^2 = ', matNorm2(elastic_L(1:3,1:3))
             WRITE(0,'(a)') 'You need to add some noise to elastic &
                &constants so as to break symmetry'
             WRITE(0,'(a)') 'Add to your input file "CVoigt_noise=1.d-4"'
            WRITE(0,'(a)') '  (do not use a noise lower than 1.d-6)'
             WRITE(0,*)
             STOP '< InitStroh >'
     ELSE IF (verbosity.GE.verbosity_max) THEN
             WRITE(out,'(a)') 'Vectors L(i,a) are linearly independant'
     END IF

     ! Check that vectors A and L are linearly independant
     matAl(1:3,1:6) = elastic_A(1:3,1:6)
     matAl(4:6,1:6) = elastic_L(1:3,1:6)
     det=MatDet(matAL)
     det_norm2 = Sqrt( Dble(det)**2 + aImag(det)**2 )
     IF ( det_norm2 .LE. 1.d-7*matNorm2(matAL)) THEN
             WRITE(0,*)
             WRITE(0,'(a)') 'WARNING < InitStroh >'
             WRITE(0,'(a)') 'Vectors A(i,a) and L(i,a) are not linearly independant'
             WRITE(0,'(a,g14.6)') '| Det(AL) |^2 = ', det_norm2
             WRITE(0,'(a,g14.6)') '| AL |^2 = ', matNorm2(matAL)
             WRITE(0,'(a)') 'You need to add some noise to elastic &
                &constants so as to break symmetry'
             WRITE(0,'(a)') 'Add to your input file "CVoigt_noise=1.d-4"'
            WRITE(0,'(a)') '  (do not use a noise lower than 1.d-6)'
             WRITE(0,*)
             STOP '< InitStroh >'
     ELSE IF (verbosity.GE.verbosity_max) THEN
             WRITE(out,'(a)') 'Vectors A(i,a) and L(i,a) are linearly independant'
     END IF

     
  END SUBROUTINE InitStroh

  !====================================================================

  SUBROUTINE Build_DStroh(burgers, force, root, elastic_A, elastic_B, elastic_L, &
        displacement, stress, angular_energy, prelog_energy_factor, interaction, out)
     ! - find vector D(n) solution of the Eq. 13.88 and 13.89 pp. 444-445
     ! - build matrices needed to calculate displacement and stress fields

     USE babel_data, ONLY : verbosity
     USE elasticity_ani
     IMPLICIT NONE
     ! Burgers vector in dislo axes: burgers(3) should be the screw component
     REAL(kind(0.d0)), dimension(1:3), intent(in) :: burgers    
     ! Force in line-force axes (it should be zero for a dislocation)
     REAL(kind(0.d0)), dimension(1:3), intent(in) :: force

     ! Outputs of subroutine Build_Elastic_A_B
     COMPLEX(kind(0.d0)), dimension(1:6), intent(in) :: root
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(in) :: elastic_A
     COMPLEX(kind(0.d0)), dimension(1:3,1:3,1:3,1:6), intent(in) :: elastic_B
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(in) :: elastic_L

     ! Matrices used to calculate displacement and stress field 
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(out) :: displacement
     COMPLEX(kind(0.d0)), dimension(1:6,1:6), intent(out) :: stress

     ! Angular dependence of energy
     REAL(kind(0.d0)), intent(out) :: angular_energy

     ! Prelogarithmic energy factor used to calculate dislocation elastic energy
     REAL(kind(0.d0)), intent(out) :: prelog_energy_factor

     ! Matrix used to calculate interaction energy between dislocations
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(out) :: interaction

     INTEGER, intent(in), optional :: out

     ! Matrix D(n) solution of the Eq. 13.88 and 13.89 pp. 444-445
     COMPLEX(kind(0.d0)), dimension(1:6) :: elastic_D

     REAL(kind(0.d0)), parameter :: pi=3.14159265358979d0
     COMPLEX(kind(0.d0)), parameter :: cmplx_i=Cmplx(0.d0,1.d0,kind(1.d0))
     COMPLEX(kind(0.d0)) :: factor, cmplx_prelog_energy_factor, Ks, pn, pm, Jp, Jpq, cmplx_E0, AL, log_factor
     INTEGER :: i, n, m
     COMPLEX(kind(0.d0)), dimension(1:3) :: u, temp
     COMPLEX(kind(0.d0)), dimension(1:6) :: sigma

     DO n=1, 6
        elastic_D(n) = -2.d0*Sum( elastic_L(1:3,n)*burgers(1:3) + elastic_A(1:3,n)*force(1:3))
     END DO

     ! Matrix used to calculate displacement and stress fields 
     !  (invert sign so as to be coherent with isotropic calculations)
     factor=-1.d0/cmplx(0.d0,2.d0*pi,kind(0.d0))
     DO n=1,6
        IF (n.EQ.4) factor=-factor
        u(1:3) = factor*elastic_A(1:3,n)*elastic_D(n)
        sigma(1) = Sum( elastic_B(1,1,1:3,n)*u(1:3) )
        sigma(2) = Sum( elastic_B(2,2,1:3,n)*u(1:3) )
        sigma(3) = Sum( elastic_B(3,3,1:3,n)*u(1:3) )
        sigma(4) = 0.5d0*Sum( (elastic_B(2,3,1:3,n)+elastic_B(3,2,1:3,n))*u(1:3) )
        sigma(5) = 0.5d0*Sum( (elastic_B(1,3,1:3,n)+elastic_B(3,1,1:3,n))*u(1:3) )
        sigma(6) = 0.5d0*Sum( (elastic_B(1,2,1:3,n)+elastic_B(2,1,1:3,n))*u(1:3) )
        displacement(1:3,n) = u(1:3)
        stress(1:6,n) = sigma(1:6)
        interaction(1:3,n) = factor*elastic_L(1:3,n)*elastic_D(n)
     END DO

     ! Prelogarithmic energy factor used to calculate dislocation elastic energy
     factor=1.d0/cmplx(0.d0,4.d0*pi,kind(0.d0))
     cmplx_prelog_energy_factor = Cmplx(0.d0, 0.d0, kind(0.d0))
     DO n=1, 6
        IF (n.EQ.4) factor=-factor
        cmplx_prelog_energy_factor = cmplx_prelog_energy_factor &
                + factor*( Sum(elastic_L(1:3,n)*burgers(1:3))**2 - Sum(elastic_A(1:3,n)*force(1:3))**2 )
     END DO

     prelog_energy_factor = Dble( cmplx_prelog_energy_factor )


     cmplx_E0 = 0.d0
     DO n=1,6
        pn=root(n)
        log_factor = 0.d0
        DO m=1,6
           AL = elastic_D(n)*Sum(elastic_A(1:3,n)*elastic_L(1:3,m) &
                - elastic_L(1:3,n)*elastic_A(1:3,m) )*elastic_D(m)
           IF (m.LE.3) THEN
                   log_factor = log_factor + AL
           ELSE
                   log_factor = log_factor - AL
                   IF (n.LE.3) THEN
                           pm=root(m)
                           cmplx_E0 = cmplx_E0 &
                                + AL/cmplx(0.d0,8.d0*pi,kind(0.d0))*log(pn-pm)
                   END IF
           END IF
        END DO
        IF (n.LE.3) THEN
                cmplx_E0 = cmplx_E0 + log_factor*log(cmplx_i+pn)/cmplx(0.d0,16.d0*pi,kind(0.d0))
        ELSE
                cmplx_E0 = cmplx_E0 + log_factor*log(cmplx_i-pn)/cmplx(0.d0,16.d0*pi,kind(0.d0))
        END IF
     END DO
     angular_energy = Dble( cmplx_E0 )

     IF ( (verbosity.GE.verbosity_max).AND.(Present(out)) ) THEN
             WRITE(out,*)
             WRITE(out,'(a)') 'Initialization for line defect using Stroh formalism'
             WRITE(out,'(a,3(1x,f0.4))') '  Burgers vector in line axes: ', &
                burgers(1:3)
             WRITE(out,'(a,3(1x,f0.4))') '  Force vector in line axes: ', &
                force(1:3)
             WRITE(out,'(a,6(g11.4,a))') '  D vector =  ( ', &
                     Dble(elastic_D(1)), ' + i*', Aimag(elastic_D(1)), ' , ', &
                     Dble(elastic_D(2)), ' + i*', Aimag(elastic_D(2)), ' , ', &
                     Dble(elastic_D(3)), ' + i*', Aimag(elastic_D(3)), ' )'
             temp(1:3) = 0.5d0*MatMul(elastic_A(:,:),elastic_D(:))
             DO i=1, 3
                WRITE(out,'(a,i1,2(a,g14.6))') &
                        '  1/2 Sum _n=1 ^6 { A_', i, '(n) D(n) } =', &
                        Dble(temp(i)), ' + i *', aImag(temp(i))
             END DO
             temp(1:3) = 0.5d0*MatMul(elastic_L(:,:),elastic_D(:))
             DO i=1, 3
                WRITE(out,'(a,i1,2(a,g14.6))') &
                        '  1/2 Sum_n=1 ^6 { L_', i, '(n) D(n) } =', &
                        Dble(temp(i)), ' + i *', aImag(temp(i))
             END DO
             WRITE(out,*)
             WRITE(out,'(a)') 'Matrix used to calculate displacement field'
             DO i=1,3
                WRITE(out,'(a,3(g11.4,a,g11.4,2x),a)') &
                     '    | ', Dble(displacement(i,1)) , ' + i*', Aimag(displacement(i,1)), &
                     Dble(displacement(i,2)) , ' + i*', Aimag(displacement(i,2)), &
                     Dble(displacement(i,3)) , ' + i*', Aimag(displacement(i,3)), ' |'
             END DO
             WRITE(out,*)
             WRITE(out,'(2(a,g13.6))') 'Angular dependence of the energy: ', &
                dble(cmplx_E0), ' + i ', aImag(cmplx_E0)
             WRITE(out,*)
             WRITE(out,'(2(a,g13.6))') 'Prelog energy factor: ', &
                dble(cmplx_prelog_energy_factor), ' + i ', aImag(cmplx_prelog_energy_factor)
             WRITE(out,*)
             DO n=1, 6
                Ks = Sum( interaction(1:3,n)*burgers(1:3) )
                WRITE(out,'(2(a,i0),2(a,g14.6))') ' -/+ b_i L_i^',n,' L_j^',n,' b_j / (2 pi i) =  ', &
                        dble(Ks), ' + i ', aImag(Ks)
             END DO
             WRITE(out,*)
             WRITE(out,*)
     END IF

  END SUBROUTINE Build_DStroh

  !====================================================================

  FUNCTION Build_inter_disloDipole(elastic_Cp, root, elastic_A, elastic_L) RESULT(Kijkl)
     ! Build matrix giving interaction energy of the core field created by
     ! dipole of dislocations with an external strain

     IMPLICIT NONE 
    ! Elastic constants in dislo or line-force axes
     REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3), intent(in) :: elastic_Cp

     ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
     !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
     !   and the equation det[a(i;j)]=0 has 6 solutions for the variable p: root(1:6) 
     COMPLEX(kind(0.d0)), dimension(1:6), intent(in) :: root
     !   corresponding vector solution of the equation a(i,j)*elastic_A(j) = 0.
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(in) :: elastic_A
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(in) :: elastic_L

     ! Interaction matrix
     REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: Kijkl


     ! Local variables and parameters
     INTEGER :: i, j, k, l, n
     COMPLEX(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: lambda
     COMPLEX :: Ix, Iy

     ! Initialization
     lambda(:,:,:,:) = Cmplx(0.d0,0.d0,kind(0.d0))

     loopt_roots: DO n=1, 6
       ! Surface integrals
       IF (n.LE.3) THEN ! Roots with positive imaginary parts
               Ix = Cmplx(0.d0,1.d0,kind(0.d0))/( Cmplx(1.d0,0.d0,kind(0.d0)) - Cmplx(0.d0,1.d0,kind(0.d0))*root(n) )
               Iy = Cmplx(1.d0,0.d0,kind(0.d0))/( Cmplx(1.d0,0.d0,kind(0.d0)) - Cmplx(0.d0,1.d0,kind(0.d0))*root(n) )
       ELSE             ! Roots with negative imaginary parts
               Ix = Cmplx(0.d0,-1.d0,kind(0.d0))/( Cmplx(1.d0,0.d0,kind(0.d0)) + Cmplx(0.d0,1.d0,kind(0.d0))*root(n) )
               Iy = Cmplx(1.d0,0.d0,kind(0.d0))/( Cmplx(1.d0,0.d0,kind(0.d0)) + Cmplx(0.d0,1.d0,kind(0.d0))*root(n) )
       END IF
       DO i=1,3; DO k=1,3
          lambda(i,1,k,1) = lambda(i,1,k,1) + elastic_A(i,n)*elastic_L(k,n)*Ix
          lambda(i,1,k,2) = lambda(i,1,k,2) + elastic_A(i,n)*elastic_L(k,n)*Ix*root(n)
          lambda(i,2,k,1) = lambda(i,2,k,1) + elastic_A(i,n)*elastic_L(k,n)*Iy
          lambda(i,2,k,2) = lambda(i,2,k,2) + elastic_A(i,n)*elastic_L(k,n)*Iy*root(n)
       END DO; END DO

     END DO loopt_roots

     ! Check that the matrix lambda is a real matrix
     DO i=1,3; DO j=1,3; DO k=1,3; DO l=1,3
        IF (Abs(aImag( lambda(i,j,k,l) )).GT. 1.d-6) THEN
                WRITE(0,'(a,4(i0,a))') 'Component lambda(',i,',',j,',',k,',',l,') is not null'
                STOP '< Build_inter_disloDipole >'
        END IF
     END DO; END DO; END DO; END DO

     !Interaction matrix with external strain
     Kijkl(:,:,:,:)=0.d0
     DO i=1,3 ; DO j=1,3 ; DO  k=1,3 ; DO l=1,3
        Kijkl(i,j,1:3,1:3) = Kijkl(i,j,1:3,1:3) + elastic_Cp(i,j,k,l) * Dble( lambda(k,l,1:3,1:3) )
     END DO ; END DO ; END DO ; END DO

     ! Interaction matrix with external stress
     Kijkl(:,:,:,:) = Dble( lambda(:,:,:,:) )   ! DEBUG

  END FUNCTION Build_inter_disloDipole
  

  !====================================================================

  FUNCTION Build_inter_lineCouple(elastic_Cp, root, elastic_A) RESULT(Kijkl)
     ! Build matrix giving interaction energy of the core field created by
     ! dipole of line-forces with an external strain

     IMPLICIT NONE 
    ! Elastic constants in dislo or line-force axes
     REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3), intent(in) :: elastic_Cp

     ! The a(i,j) matrix of Eq. 13.57 p. 438 is given by
     !   a(i,j) = a_poly(0,i,j) + a_poly(1,i,j)*p + a_poly(2,i,j)*p*p
     !   and the equation det[a(i;j)]=0 has 6 solutions for the variable p: root(1:6) 
     COMPLEX(kind(0.d0)), dimension(1:6), intent(in) :: root
     !   corresponding vector solution of the equation a(i,j)*elastic_A(j) = 0.
     COMPLEX(kind(0.d0)), dimension(1:3,1:6), intent(in) :: elastic_A

     ! Interaction matrix
     REAL(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: Kijkl


     ! Local variables and parameters
     INTEGER :: i, j, k, l, n
     COMPLEX(kind(0.d0)), dimension(1:3,1:3,1:3,1:3) :: lambda
     COMPLEX :: Ix, Iy

     !Interaction matrix with external strain
     Kijkl(:,:,:,:)=0.d0
     DO i=1,3 ; DO j=1,3 ; DO  k=1,3 ; DO l=1,3
        Kijkl(i,j,i,j) = 1.d0
     END DO ; END DO ; END DO ; END DO
     RETURN

     !==========================
     ! Initialization
     lambda(:,:,:,:) = Cmplx(0.d0,0.d0,kind(0.d0))

     loopt_roots: DO n=1, 6
       ! Surface integrals
       IF (n.LE.3) THEN ! Roots with positive imaginary parts
               Ix = Cmplx(0.d0,1.d0,kind(0.d0))/( Cmplx(1.d0,0.d0,kind(0.d0)) - Cmplx(0.d0,1.d0,kind(0.d0))*root(n) )
               Iy = Cmplx(1.d0,0.d0,kind(0.d0))/( Cmplx(1.d0,0.d0,kind(0.d0)) - Cmplx(0.d0,1.d0,kind(0.d0))*root(n) )
       ELSE             ! Roots with negative imaginary parts
               Ix = Cmplx(0.d0,-1.d0,kind(0.d0))/( Cmplx(1.d0,0.d0,kind(0.d0)) + Cmplx(0.d0,1.d0,kind(0.d0))*root(n) )
               Iy = Cmplx(1.d0,0.d0,kind(0.d0))/( Cmplx(1.d0,0.d0,kind(0.d0)) + Cmplx(0.d0,1.d0,kind(0.d0))*root(n) )
       END IF
       DO i=1,3; DO k=1,3
          lambda(i,1,k,1) = lambda(i,1,k,1) + elastic_A(i,n)*elastic_A(k,n)*Ix
          lambda(i,1,k,2) = lambda(i,1,k,2) + elastic_A(i,n)*elastic_A(k,n)*Ix*root(n)
          lambda(i,2,k,1) = lambda(i,2,k,1) + elastic_A(i,n)*elastic_A(k,n)*Iy
          lambda(i,2,k,2) = lambda(i,2,k,2) + elastic_A(i,n)*elastic_A(k,n)*Iy*root(n)
       END DO; END DO

     END DO loopt_roots

     ! Check that the matrix lambda is a real matrix
     DO i=1,3; DO j=1,3; DO k=1,3; DO l=1,3
        IF (Abs(aImag( lambda(i,j,k,l) )).GT. 1.d-6) THEN
                WRITE(0,'(a,4(i0,a))') 'Component lambda(',i,',',j,',',k,',',l,') is not null'
                STOP '< Build_inter_lineCouple >'
        END IF
     END DO; END DO; END DO; END DO

     ! Interaction matrix with external strain
     Kijkl(:,:,:,:)=0.d0
     DO i=1,3 ; DO j=1,3 ; DO  k=1,3 ; DO l=1,3
        Kijkl(i,j,1:3,1:3) = Kijkl(i,j,1:3,1:3) + elastic_Cp(i,j,k,l) * Dble( lambda(k,l,1:3,1:3) )
     END DO ; END DO ; END DO ; END DO

     ! Interaction matrix with external stress
     !!$Kijkl(:,:,:,:) = Dble( lambda(:,:,:,:) )   ! DEBUG

  END FUNCTION Build_inter_lineCouple
  

  !====================================================================

END MODULE elasticity_Stroh
