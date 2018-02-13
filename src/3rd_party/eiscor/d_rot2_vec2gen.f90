#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2_vec2gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generator for a Givens rotation represented by 2
! real numbers: A strictly real cosine and a scrictly real sine.  The first
! column is constructed to be parallel with the real vector [A,B]^T. 
!
! The rot2 refers to a rotation desribed by two double and the vec2 to a vector
! of lenght two described by two double.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  A               REAL(8) 
!                    the first component of the real vector [A,B]^T
!
!  B               REAL(8) 
!                    the second component of the real vector [A,B]^T
!
! OUTPUT VARIABLES:
!
!  C               REAL(8)
!                    on exit contains the generator for the cosine
!                    component of the Givens rotation
!
!  S               REAL(8)
!                    on exit contains the generator for the sine
!                    component of the Givens rotation
!
!  NRM             REAL(8)
!                    on exit contains the norm of the vector [A,B]^T
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_vec2gen(A,B,C,S,NRM)

  implicit none
  
  ! input variables
  real(8), intent(in) :: A,B
  real(8), intent(inout) :: C,S,NRM

  ! construct rotation
  if ((A.EQ.0d0).AND.(B.EQ.0d0)) then
     C = 1d0
     S = 0d0
     NRM = 0d0
  else if (abs(A).GE.abs(B)) then
     S = B/A
     NRM = sign(sqrt(1d0 + S*S),A)
     C =  1.d0/NRM
     S =  S*C
     NRM =  A*NRM
  else
     C = A/B;
     NRM = sign(sqrt(1d0 + C*C),B)
     S =  1.d0/NRM
     C =  C*S
     NRM =  B*NRM
  end if
           
end subroutine d_rot2_vec2gen
