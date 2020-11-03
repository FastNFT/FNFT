#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_mergebulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the product of two real Given's rotations and 
! and stores the output in one of the input arrays.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  TOP             LOGICAL
!                    .TRUE.: merge bulge at top
!                    .FALSE.: merge bulge at bottom
!
!  Q               REAL(8) array of dimension (2)
!                    generator unitary matrix
!
!  B               REAL(8) array of dimension (2)
!                    generator for bulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_mergebulge(TOP,Q,B)

  implicit none

  ! input variables
  logical, intent(in) :: TOP
  real(8), intent(inout) :: Q(2), B(2)

  ! compute variables
  real(8) :: c, s, nrm

  ! merge at top
  if(TOP)then

    ! compute new generators
    c = B(1)*Q(1)-B(2)*Q(2)
    s = B(2)*Q(1)+B(1)*Q(2)
    
    ! compute new generators
    call d_rot2_vec2gen(c,s,Q(1),Q(2),nrm)
    
  ! merge at bottom
  else

    ! compute new generators
    c = Q(1)*B(1)-Q(2)*B(2)
    s = Q(2)*B(1)+Q(1)*B(2)
    
    ! compute new generators
    call d_rot2_vec2gen(c,s,Q(1),Q(2),nrm)
    
  end if

end subroutine d_orthfact_mergebulge
