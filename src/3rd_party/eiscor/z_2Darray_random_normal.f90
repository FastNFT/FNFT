#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_2Darray_random_normal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine provides a matrix of pseudo-random complex numbers independent, 
! standard, normally distributed (expected value 0, standard deviation 1).
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
!  A               COMPLEX(8) array of dimension (M,N)
!                    matrix of complex random numbers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_2Darray_random_normal(M,N,A)

  implicit none
  
  ! input variables
  integer, intent(in) :: M, N
  complex(8), intent(inout) :: A(M,N)
  
  ! compute variables
  integer :: ii, jj
  real(8) :: u, v
  
  ! check N
  if ((N .LT. 1).OR.(M .LT. 1)) then
    return
  end if
  
  do jj = 1,M
     do ii = 1,N
        
        call z_scalar_random_normal(u,v)
        A(jj,ii) = cmplx(u,v,kind=8)
        
     end do
  end do
  
end subroutine z_2Darray_random_normal
