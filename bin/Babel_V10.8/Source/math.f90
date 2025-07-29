MODULE Math

  INTERFACE matDet
          ! Calculate determinant of a real or complex n*n matrix
          MODULE PROCEDURE real_matDet, complex_matDet
  END INTERFACE matDet

  INTERFACE vecNorm2
          ! Calculate Norm 2 of a real or complex n vector
          MODULE PROCEDURE real_vecNorm2, complex_vecNorm2
  END INTERFACE vecNorm2

  INTERFACE matNorm2
          ! Calculate Norm 2 of a real or complex n*n matrix
          MODULE PROCEDURE real_matNorm2, complex_matNorm2
  END INTERFACE matNorm2

CONTAINS

  FUNCTION CrossProduct(ex,ey) RESULT(ez)
    
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: ex, ey
    REAL(kind(0.d0)), dimension(1:3) :: ez

    ez(1) = ex(2)*ey(3) - ex(3)*ey(2)
    ez(2) = ex(3)*ey(1) - ex(1)*ey(3)
    ez(3) = ex(1)*ey(2) - ex(2)*ey(1)

  END FUNCTION CrossProduct

  FUNCTION Colinear(ex,ey) RESULT(test)
    
    USE Babel_Data
    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: ex, ey
    LOGICAL :: test

    REAL(kind(0.d0)) :: n2

    n2 = ( ex(2)*ey(3) - ex(3)*ey(2) )**2 &
        + ( ex(3)*ey(1) - ex(1)*ey(3) )**2 &
        + ( ex(1)*ey(2) - ex(2)*ey(1) )**2
    test = (n2.LE.distance_Zero2)

  END FUNCTION Colinear

  FUNCTION ScalarTripleProduct(ex,ey,ez) RESULT(stp)

    IMPLICIT NONE
    REAL(kind(0.d0)), dimension(1:3), intent(in) :: ex, ey, ez
    REAL(kind(0.d0)) :: stp

    stp = ez(1) *( ex(2)*ey(3) - ex(3)*ey(2) ) &
        + ez(2) *( ex(3)*ey(1) - ex(1)*ey(3) ) &
        + ez(3) *( ex(1)*ey(2) - ex(2)*ey(1) )

  END FUNCTION ScalarTripleProduct

  FUNCTION real_vecNorm2(A) result(n2)

    implicit none

    REAL(kind(0.d0)), dimension(:), intent(in) :: A
    REAL(kind(0.d0)) :: n2

    n2 = Sqrt( Sum( A(:)**2 ) )

  END FUNCTION real_vecNorm2

  FUNCTION complex_vecNorm2(A) result(n2)

    implicit none

    COMPLEX(kind(0.d0)), dimension(:), intent(in) :: A
    REAL(kind(0.d0)) :: n2

    n2 = Sqrt( Sum( Dble( A(:) )**2 + aImag( A(:) )**2 ) )

  END FUNCTION complex_vecNorm2

  FUNCTION real_matNorm2(A) result(n2)

    implicit none

    REAL(kind(0.d0)), dimension(:,:), intent(in) :: A
    REAL(kind(0.d0)) :: n2

    n2 = Sqrt( Sum( A(:,:)**2 ) )

  END FUNCTION real_matNorm2

  FUNCTION complex_matNorm2(A) result(n2)

    implicit none

    COMPLEX(kind(0.d0)), dimension(:,:), intent(in) :: A
    REAL(kind(0.d0)) :: n2

    n2 = Sqrt( Sum( Dble( A(:,:) )**2 + aImag( A(:,:) )**2 ) )

  END FUNCTION complex_matNorm2

  FUNCTION iNormInf(z) RESULT(norm)

    IMPLICIT NONE

    COMPLEX(kind(0.d0)), intent(in) :: z
    REAL(kind(0.d0)) :: norm

    norm = Abs( dble(z) ) + Abs( aImag(z) )

  END FUNCTION iNormInf

  FUNCTION real_matDet(A) result(det)
    ! Determinant of a real matrix

    implicit none

    REAL(kind(0.d0)), dimension(:,:), intent(in) :: A
    REAL(kind(0.d0)) :: det
    INTEGER :: nl, nc, n, info, nPermut
    REAL(kind(0.d0)), dimension(:,:), allocatable :: Atemp
    INTEGER, dimension(:), allocatable :: iPiv
    EXTERNAL :: DGETRF

    ! Number of lines and columns
    nl = size(A,1) ; nc=size(A,2)
    IF ( nl.NE.nc ) THEN
            WRITE(0,'(a,i0)') 'Number of lines:   nl=', nl
            WRITE(0,'(a,i0)') 'Number of columns: nc=', nc
            STOP '< real_matDet >: matrix is not square'
    END IF

    SELECT CASE(nl)
    CASE(1)
            det = a(1,1)
    CASE(2)
            det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    CASE(3)
            det =  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
                 + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1) &
                 - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3)
    CASE DEFAULT
            ! Use Lapack subroutine to make LU factorization
            ALLOCATE(Atemp(nl,nl))
            ALLOCATE(iPiv(nl))
            Atemp(:,:) = A(:,:)
            CALL DGETRF(nl, nl, Atemp, nl, iPiv, info)
            IF (info.LT.0) THEN
                    WRITE(0,'(a,i0)') 'Info = ', info
                    WRITE(0,'(a)') 'Error when calling Lapack subroutine DGETRF'
                    WRITE(0,'(a)') 'if INFO = -i, the i-th argument had an illegal value'
                    STOP '< Real_matDet >:'
            ELSE IF (info.GT.0) THEN
                    det=0.d0
            ELSE
                    ! Determinant of the U matrix
                    det = 1.d0
                    DO n=1, nl
                       det = det*Atemp(n,n)
                    END DO
                    ! Number of permutations
                    nPermut=0
                    DO n=1, nl
                       IF (iPiv(n).NE.n) nPermut=nPermut+1
                    END DO
                    IF (Modulo(nPermut,2).EQ.1) det=-det
            END IF
            DEALLOCATE(Atemp, iPiv)
    END SELECT


  END FUNCTION real_matDet

  FUNCTION complex_matDet(A) result(det)

    implicit none

    COMPLEX(kind(0.d0)), dimension(:,:), intent(in) :: A
    COMPLEX(kind(0.d0)) :: det
    INTEGER :: nl, nc, n, info, nPermut
    COMPLEX(kind(0.d0)), dimension(:,:), allocatable :: Atemp
    INTEGER, dimension(:), allocatable :: iPiv
    EXTERNAL :: ZGETRF


    ! Number of lines and columns
    nl = size(A,1) ; nc=size(A,2)
    IF ( nl.NE.nc ) THEN
            WRITE(0,'(a,i0)') 'Number of lines:   nl=', nl
            WRITE(0,'(a,i0)') 'Number of columns: nc=', nc
            STOP '< complex_matDet >: matrix is not square'
    END IF

    !!$WRITE(0,'(a)') 'Call to complex_matDet'             ! DEBUG
    !!$WRITE(0,'(2(a,i0))') '  nl=', nl, ' nc=', nc        ! DEBUG

    SELECT CASE(nl)
    CASE(1)
            det = a(1,1)
    CASE(2)
            det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    CASE(3)
            det =  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
                 + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1) &
                 - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3)
    CASE DEFAULT
            ! Use Lapack subroutine to make LU factorization
            ALLOCATE(Atemp(nl,nl))
            ALLOCATE(iPiv(nl))
            Atemp(:,:) = A(:,:)
            CALL ZGETRF(nl, nl, Atemp, nl, iPiv, info)
            IF (info.LT.0) THEN
                    WRITE(0,'(a,i0)') 'Info = ', info
                    WRITE(0,'(a)') 'Error when calling Lapack subroutine ZGETRF'
                    WRITE(0,'(a)') 'if INFO = -i, the i-th argument had an illegal value'
                    STOP '< Real_matDet >:'
            ELSE IF (info.GT.0) THEN
                    det = Cmplx(0.d0, 0.d0, kind(0.d0))
            ELSE
                    ! Determinant of the U matrix
                    det = Cmplx(1.d0, 0.d0, kind(0.d0))
                    DO n=1, nl
                       det = det*Atemp(n,n)
                       !!$WRITE(0,'(2(a,i0),2(a,g14.6))') ' U(', n,',',n,') = ', &     ! DEBUG
                                !!$Dble(Atemp(n,n)), ' + i * ', aImag(Atemp(n,n))      ! DEBUG
                    END DO
                    ! Number of permutations
                    nPermut=0
                    DO n=1, nl
                       IF (iPiv(n).NE.n) nPermut=nPermut+1
                    END DO
                    IF (Modulo(nPermut,2).EQ.1) det=-det
            END IF
            DEALLOCATE(Atemp, iPiv)
    END SELECT


  END FUNCTION complex_matDet

  SUBROUTINE matinv(A, B)
    ! return matrix B which is the inverse of matrix A

    implicit none

    REAL(kind(0.d0)), dimension(:,:), intent(in) :: A
    REAL(kind(0.d0)), dimension(:,:), intent(out) :: B

    REAL(kind(0.d0)) :: invdet
    EXTERNAL :: DGESV

    ! Variables used with Lapack subroutine
    INTEGER :: i, n, INFO
    INTEGER, dimension(:), allocatable :: IPIV
    REAL(kind(0.d0)), dimension(:,:), allocatable :: tempA, invA

    if (size(A,1).NE.size(A,2)) &
         STOP '< MatInv >: matrix to invert is not square'

    SELECT CASE(size(A,1))
    CASE(1)
       B(1,1)=1.d0/A(1,1)
       
    CASE(2)
       invdet=1.d0/( A(1,1)*A(2,2)-A(1,2)*A(2,1) )
       B(1,1)=A(2,2)*invdet
       B(2,2)=A(1,1)*invdet
       B(1,2)=-A(1,2)*invdet
       B(2,1)=-A(2,1)*invdet
       
    CASE(3)
       invdet=1.d0/real_matDet(A)
       
       b(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
       b(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
       b(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
       
       b(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
       b(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
       b(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)
       
       b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
       b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
       b(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

       b(:,:)=b(:,:)*invdet
       
    CASE DEFAULT
       ! Use Lapack subroutine DGESV       
      n=size(A,1)
      allocate(tempA(n,n))
      tempA(1:n,1:n)=A(1:n,1:n) ! Needed in order to not modify A
      allocate(IPIV(1:n))
      allocate(invA(n,n))
      invA(1:n,1:n)=0.d0
      do i=1,n
         invA(i,i)=1.d0
      enddo
      call DGESV(n, n, tempA, n, IPIV, invA, n, INFO)
      if (INFO.NE.0) then
         write(0,*) 'Result of DGESV subroutine: INFO=', INFO 
         write(0,*)' < 0: if INFO = -i, the i-th argument had an illegal value'
         write(0,*)' > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization'
         write(0,*)'has been completed, but the factor U is exactly'
         write(0,*)'singular, so the solution could not be computed.'
         STOP '< MatInv >'
      endif
      B(1:n,1:n)=invA(1:n,1:n)
      deallocate(tempA, invA, IPIV)

    END SELECT

  END SUBROUTINE matinv


END MODULE Math
