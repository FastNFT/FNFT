#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2_swapdiag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes the generator for a Givens rotation represented
! by 2 real numbers: a strictly real cosine and a scrictly real sine 
! through a 2x2 diagonal matrix with entries +/-1.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  D               REAL(8) array of dimension (2)
!                    array of generators for complex diagonal matrix
!                    the entries of D must be +/-1
!
!  B               REAL(8) array of dimension (2)
!                    generator for a Givens rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_swapdiag(D,B)

  implicit none
  
  ! input variables
  real(8), intent(in) :: D(2)
  real(8), intent(inout) :: B(2)
  
  ! change sign in D not scalar
  if (D(1).NE.D(2)) then
    B(2) = -B(2)
  end if 
  
end subroutine d_rot2_swapdiag
