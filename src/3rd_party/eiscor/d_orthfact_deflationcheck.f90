#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_deflationcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in an orthogonal upper hessenberg 
! matrix that is stored as a product of givens rotations and a complex 
! diagonal matrix. When a deflation occurs the corresponding rotation
! is set to the identity matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (N)
!                    array of generators for complex diagonal matrix
!                    entries must be +/-1
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     index of the last known deflation
!                     on output contains index of newest deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_deflationcheck(N,Q,D,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(2*(N-1)), D(N)
  integer, intent(inout) :: ZERO

  ! compute variables
  integer :: ii
  real(8), parameter :: tol = EISCOR_DBL_EPS

  ! initialize ZERO
  ZERO = 0
  
  ! check for deflation
  do ii=1,(N-1)
  
    ! deflate if subdiagonal is small enough
    if(abs(Q(2*(N-ii))) < tol)then
     
      ! update D
      D(N-ii) = sign(1d0,Q(2*(N-ii)-1)*D(N-ii))
      D(N-ii+1) = sign(1d0,Q(2*(N-ii)-1)*D(N-ii+1))
 
      ! update subsequent Q
      if (ii > 1) then
        Q(2*(N-ii+1)) = sign(1d0,Q(2*(N-ii)-1))*Q(2*(N-ii+1))
      end if

      ! set Q to identity
      Q(2*(N-ii)-1) = 1d0
      Q(2*(N-ii)) = 0d0
       
      ! set ZERO
      ZERO = max(0,N-ii)
        
      ! exit loop  
      exit

    end if

  end do
  
end subroutine d_orthfact_deflationcheck
