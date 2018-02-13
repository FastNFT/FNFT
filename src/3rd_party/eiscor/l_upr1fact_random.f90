#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! l_upr1fact_random
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This function returns the new position flag for a single shift 
! iteration on a twisted upr1 pencil based on the position flags from
! the current factorization. The returned flag is chosen randomly. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function l_upr1fact_random(N,P)

  implicit none
  
  ! return type
  logical :: l_upr1fact_random
  
  ! input variables
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)

  ! compute variables
  real(8) :: num

  ! call random number
  call random_number(num)
  
  ! set output
  if (num < 5d-1) then
    l_upr1fact_random = .FALSE.
  else
    l_upr1fact_random = .TRUE.
  end if

end function l_upr1fact_random
