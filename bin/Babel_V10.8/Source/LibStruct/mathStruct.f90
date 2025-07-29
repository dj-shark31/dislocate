MODULE MathStruct

  IMPLICIT NONE

CONTAINS

  FUNCTION mat3Norm2(A) result(n2)

    implicit none

    REAL(kind(0.d0)), dimension(3,3), intent(in) :: A
    REAL(kind(0.d0)) :: n2

    n2 = Sqrt( Sum( A(:,:)**2 ) )

  END FUNCTION mat3Norm2

  SUBROUTINE mat3inv(A, B)
    ! return matrix B which is the inverse of matrix A
    ! A and B are 3x3 matrices

    implicit none

    REAL(kind(0.d0)), dimension(3,3), intent(in) :: A
    REAL(kind(0.d0)), dimension(3,3), intent(out) :: B

    REAL(kind(0.d0)) :: invdet

       
      b(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
      b(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
      b(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
      
      b(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
      b(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
      b(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)
      
      b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
      b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
      b(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

      invdet = 1.d0/( a(1,1)*b(1,1) + a(1,2)*b(2,1) + a(1,3)*b(3,1) )
      b(:,:)=b(:,:)*invdet

  END SUBROUTINE mat3inv

  FUNCTION mat3Det(A) result(det)
    ! Determinant of a real 3x3 matrix

    implicit none

    REAL(kind(0.d0)), dimension(3,3), intent(in) :: A
    REAL(kind(0.d0)) :: det

            det =  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
                 + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1) &
                 - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3)

  END FUNCTION mat3Det

END MODULE MathStruct

