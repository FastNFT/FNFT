#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_unitvec3gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generators for a complex rotation represented
! by 3 real numbers: the real and imaginary parts of a complex cosine,
! CR and CI, and a strictly real sine, S. The CR, CI and S are 
! constructed to be parallel with the vector [AR+iAI,B]^T, where AR, 
! AI and B are real and i = sqrt(-1).
!
! The vector [AR+iAI,B]^T is assumed to be nearly unit, that is 
!
!   abs(AR^2 + AI^2 + B^2 - 1) = O(eps).
! 
! In this case it is possible to compute the norm without a squareroot using a
! single newton iteration.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  AR, AI          REAL(8) 
!                    real and imaginary part of the first component 
!                    of the complex vector [AR+iAI,B]^T
!
!  B               REAL(8) 
!                    the second component 
!                    of the complex vector [AR+iAI,B]^T
!
! OUTPUT VARIABLES:
!
!  CR, CI          REAL(8)
!                    on exit contains the generators for the cosine
!                    component of the Givens rotation
!
!  S               REAL(8)
!                    on exit contains the generator for the sine
!                    component of the Givens rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_unitvec3gen(AR,AI,B,CR,CI,S)

  implicit none
  
  ! input variables
  real(8), intent(in) :: AR, AI, B
  real(8), intent(inout) :: CR, CI, S
  
  ! compute variables
  real(8) :: X

  ! compute difference from unity and approximate 1/sqrt(x)
  X = 5d-1*(3d0 - (AR*AR + AI*AI + B*B))
  
  ! construct rotation
  CR = AR*X
  CI = AI*X
  S  =  B*X
    
end subroutine z_rot3_unitvec3gen
