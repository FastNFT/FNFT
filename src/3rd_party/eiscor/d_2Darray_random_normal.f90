#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_2Darray_random_normal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine provides a matrix of pseudo-random number independent, standard,
! normally distributed (expected value 0, standard deviation 1).
! Call u_randomseed_initialize(INFO) first to set the seed of the number 
! generator.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  M, N            INTEGER
!                    positive integers, dimension of A
!
!
! OUTPUT VARIABLES:
!
!  A               REAL(8) array of dimension (M,N)
!                    matrix of real random numbers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_2Darray_random_normal(M,N,A)

  implicit none
  
  ! input variables
  integer, intent(in) :: M, N
  real(8), intent(inout) :: A(M,N)
  
  ! compute variables
  integer :: ii, jj
  
  ! check N
  if ((N .LT. 1).OR.(M .LT. 1)) then
    return
  end if
  
  do jj = 1,M
     do ii = 1,N
        
        call d_scalar_random_normal(A(jj,ii))
        
     end do
  end do
  
end subroutine d_2Darray_random_normal
