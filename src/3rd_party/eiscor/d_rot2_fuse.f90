#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2_fuse
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
!  DIR             LOGICAL
!                    .TRUE.: A = A*B
!                    .FALSE.: B = A*B
!
!  A,B             REAL(8) array of dimension (2)
!                    generators for rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_fuse(DIR,A,B)

  implicit none

  ! input variables
  logical, intent(in) :: DIR
  real(8), intent(inout) :: A(2), B(2)

  ! compute variables
  real(8) :: c, s, nrm

  ! merge at top
  if(DIR)then

    ! compute new generators
    c = A(1)*B(1)-A(2)*B(2)
    s = A(2)*B(1)+A(1)*B(2)
    
    ! compute new generators
    call d_rot2_vec2gen(c,s,A(1),A(2),nrm)
    
  ! merge at bottom
  else

    ! compute new generators
    c = A(1)*B(1)-A(2)*B(2)
    s = A(2)*B(1)+A(1)*B(2)
    
    ! compute new generators
    call d_rot2_vec2gen(c,s,B(1),B(2),nrm)
    
  end if

end subroutine d_rot2_fuse
