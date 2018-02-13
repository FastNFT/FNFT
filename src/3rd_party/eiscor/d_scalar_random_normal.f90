#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_scalar_random_normal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine provides a pseudo-random number independent, standard,
! normally distributed (expected value 0, standard deviation 1). The routine
! generates two uniformly distributed random numbers using Fortran's built-in
! routine and then uses the Box-Muller transform. 
!
! Call u_randomseed_initialize(INFO) first to set the seed of the number 
! generator.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! OUTPUT VARIABLES:
!
!  NUM             REAL(8) 
!                    random numbers 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_scalar_random_normal(NUM)

  implicit none
  
  ! input variables
  real(8), intent(inout) :: NUM

  ! compute variables
  double precision :: u,v,s,pi = EISCOR_DBL_PI
  integer :: ii
 
  do ii = 1, 100
     call random_number(u)
     call random_number(v)
     s = u**2 + v**2
              
     ! Box-Muller transform
     if ((s > 0d0) .AND. (s < 1d0)) then
        NUM = dcos(2.d0*pi*u)*dsqrt(-2.d0*dlog(v))
        exit
     end if
  end do

end subroutine d_scalar_random_normal
