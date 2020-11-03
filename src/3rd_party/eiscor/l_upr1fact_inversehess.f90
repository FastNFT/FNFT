#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! l_upr1fact_inversehess
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This function returns the new position flag for a single shift 
! iteration on a twisted upr1 pencil based on the position flags from
! the current factorization. The returned flag is always .TRUE. 
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
function l_upr1fact_inversehess(N,P)

  implicit none
  
  ! return type
  logical :: l_upr1fact_inversehess
  
  ! input variables
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  
  ! set output
  l_upr1fact_inversehess = .TRUE.

end function l_upr1fact_inversehess
