#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_1Darray_random_normal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine provides a vector of pseudo-random complex numbers independent, 
! standard, normally distributed (expected value 0, standard deviation 1).
! Call u_randomseed_initialize(INFO) first to set the seed of the number 
! generator.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    positive integers, dimension of A
!
!
! OUTPUT VARIABLES:
!
!  A               COMPLEX(8) array of dimension (N)
!                    vector of complex random numbers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_1Darray_random_normal(N,A)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(inout) :: A(N)
  
  ! compute variables
  integer :: ii
  real(8) :: u, v
  
  ! check N
  if (N .LT. 1) then
    return
  end if
  
  do ii = 1,N
     
     call z_scalar_random_normal(u,v)
     A(ii) = cmplx(u,v,kind=8)
        
  end do
  
end subroutine z_1Darray_random_normal
