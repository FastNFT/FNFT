#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_2x2diagblock 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes either the top or bottom 2x2 diagonal block 
! of an orthogonal upper-hessenberg matrix that is stored as a product 
! of givens rotations and a diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  TOP             LOGICAL
!                    .TRUE.: top block is computed
!                    .FALSE.: bottom block is computed
!
!  Q               REAL(8) array of dimension (4)
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (2)
!                    array of generators for diagonal matrix
!                    entries must be +/-1
!
! OUTPUT VARIABLES:
!
!  H               REAL(8) array of dimension (2,2)
!                    on exit contains the desired 2x2 block
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_2x2diagblock(TOP,Q,D,H)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: TOP
  real(8), intent(in) :: Q(4), D(2)
  real(8), intent(inout) :: H(2,2)
  
  ! TOP == .TRUE.
  if (TOP) then

    ! initialize H
    H(1,1) = Q(1)*D(1)
    H(2,1) = Q(2)*D(1)
    H(1,2) = -Q(2)*Q(3)*D(2)
    H(2,2) = Q(1)*Q(3)*D(2)

  ! TOP == .FALSE.
  else

    ! initialize H
    H(1,1) = Q(3)*Q(1)*D(1)
    H(2,1) = Q(4)*D(1)
    H(1,2) = -Q(4)*Q(1)*D(2)
    H(2,2) = Q(3)*D(2)

  end if 

end subroutine d_orthfact_2x2diagblock
