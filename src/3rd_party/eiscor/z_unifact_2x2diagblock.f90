#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_2x2diagblock
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a two by two diagonal block of a unitary upper 
! hessenberg matrix that is stored as a product of givens rotations 
! and a complex diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  TOP             LOGICAL
!                    .TRUE.: top block is computed
!                    .FALSE.: bottom block is computed
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (4)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  H              COMPLEX(8) array of dimension (2,2)
!                   on exit contains the desired 2x2 block
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_2x2diagblock(TOP,Q,D,H)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: TOP
  real(8), intent(in) :: Q(6), D(4)
  complex(8), intent(inout) :: H(2,2)
 
  ! TOP == .TRUE.
  if (TOP) then
 
    ! initialize H
    H(1,1) = cmplx(Q(1),Q(2),kind=8)
    H(2,1) = cmplx(Q(3),0d0,kind=8)
    H(1,2) = -H(2,1)
    H(2,2) = conjg(H(1,1))
   
    ! second rotation 
    H(:,2) = H(:,2)*cmplx(Q(4),Q(5),kind=8)

  ! TOP == .FALSE.
  else
 
    ! initialize H
    H(1,1) = cmplx(Q(4),Q(5),kind=8)
    H(2,1) = cmplx(Q(6),0d0,kind=8)
    H(1,2) = -H(2,1)
    H(2,2) = conjg(H(1,1))
   
    ! second rotation
    H(1,:) = H(1,:)*cmplx(Q(1),-Q(2),kind=8)
    
  end if

  ! diagaonal
  H(:,1) = H(:,1)*cmplx(D(1),D(2),kind=8)
  H(:,2) = H(:,2)*cmplx(D(3),D(4),kind=8)
  
end subroutine z_unifact_2x2diagblock
