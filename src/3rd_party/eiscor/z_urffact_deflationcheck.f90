#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine diagonalizes a unitary upper hessenberg matrix that is stored as
! a product of N Givens rotations, without computing square roots.
!
! | u1       -v1 |
! | v1  conj(u1) | | u2       -v2 | 
!                  | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                   | v3  conj(u3) | | 0  u4 |                                   
!                                                                     
! The square root free algorithm only requires the storage of the vi^2,
! so the arrays U and VV contain the following:
!
!  U(i) = ui
! VV(i) = vi^2
!
! The input must satisfy the following:
!
!  |U(i)|^2 + VV(i) = 1
!             VV(N) = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for Givens rotations
!
!  VV              REAL(8) array of dimension N
!                    array of real generators (squared) for Givens rotations
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                    largest index such that VV(i) < tol
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_deflationcheck(N,U,VV,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ZERO
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)

  ! compute variables
  integer :: ii
  real(8), parameter :: tol = (EISCOR_DBL_EPS)**2
  real(8) :: xx

  ! intialize ZERO
  ZERO = 0
  
  ! check for deflation
  do ii=1,N
  
    ! deflate if subdiagonal is small enough
    if (VV(N+1-ii) < tol) then
        
      ! set ZERO
      ZERO = N+1-ii

      ! set rotation to diagonal
      VV(ZERO) = 0d0
        
      ! renormalize U
      xx = dble(U(ZERO))**2 + aimag(U(ZERO))**2
      U(ZERO) = 5d-1*U(ZERO)*(3d0-xx)

    end if

  end do
  
end subroutine z_urffact_deflationcheck
