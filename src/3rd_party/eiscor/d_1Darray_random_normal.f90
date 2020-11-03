#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_1Darray_random_normal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine provides a vector of pseudo-random number independent, standard,
! normally distributed (expected value 0, standard deviation 1).
! Call u_randomseed_initialize(INFO) first to set the seed of the number 
! generator.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    positive integer, dimension of A
!
!
! OUTPUT VARIABLES:
!
!  A               REAL(8) array of dimension (N)
!                    vector of real random numbers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_1Darray_random_normal(N,A)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(inout) :: A(N)
  
  ! compute variables
  integer :: ii
  
  ! check N
  if (N < 1) then
    return
  end if
  
  do ii = 1,N

    call d_scalar_random_normal(A(ii))
    
  end do

end subroutine d_1Darray_random_normal
