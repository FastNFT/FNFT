#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first transformation in a single shift
! iteration for a unitary upper hessenberg matix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (4)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
!  SHFT            COMPLEX(8) 
!                    contains the shift need for the first transformation
!
! OUTPUT VARIABLES:
!
!  B               REAL(8) array of dimension (3)
!                    on exit contains the generators for the first
!                    transformation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_buildbulge(Q,D,SHFT,B)
  
  implicit none
  
  ! input variables
  real(8), intent(in) :: Q(6), D(4)
  complex(8), intent(in) :: SHFT
  real(8), intent(inout) :: B(3)
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: block(2,2)
  
  ! get top block
  call z_unifact_2x2diagblock(.TRUE.,Q,D,block) 
      
  ! shift first entry
  block(1,1) = block(1,1) - SHFT
  
  ! bulge
  call z_rot3_vec4gen(dble(block(1,1)),aimag(block(1,1)),dble(block(2,1)), &
    aimag(block(2,1)),B(1),B(2),B(3),nrm)
      
end subroutine z_unifact_buildbulge
