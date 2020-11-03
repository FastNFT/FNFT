#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine uses the generators for three Givens rotations represented
! by 3 real numbers: the real and imaginary parts of a complex cosine
! and a scrictly real sine and performs a turnover.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  Q1, Q2, Q3       REAL(8) arrays of dimension (2)
!                    generators for givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_turnover(Q1,Q2,Q3)

  implicit none
  
  ! input variables
  real(8), intent(inout) :: Q1(2), Q2(2), Q3(2)

  ! compute variables
  real(8) :: nrm
  real(8) :: a, b 
  real(8) :: c1, s1
  real(8) :: c2, s2
  real(8) :: c3, s3
  real(8) :: c4, s4
  real(8) :: c5, s5
  real(8) :: c6, s6
  
  ! set local variables
  c1 = Q1(1)
  s1 = Q1(2)
  c2 = Q2(1)
  s2 = Q2(2)
  c3 = Q3(1)
  s3 = Q3(2)
  
  ! compute first rotation
  a = s1*c3 + c1*c2*s3
  b = s2*s3
  call d_rot2_vec2gen(a,b,c4,s4,nrm)

  ! compute second rotation
  a = c1*c3 - s1*c2*s3
  b = nrm
  call d_rot2_unitvec2gen(a,b,c5,s5)

  ! compute third rotation
  ! s5 == 0 
  a = c1*s2*s4 + c2*c4
  if (s5.EQ.0d0) then
    b = c5*c1*c4*s2
  ! all other cases
  else
    b = s1*s2/s5
  end if  
  call d_rot2_unitvec2gen(a,b,c6,s6)

  ! set output
  Q1(1) = c5
  Q1(2) = s5
  Q2(1) = c6
  Q2(2) = s6
  Q3(1) = c4
  Q3(2) = s4
       
end subroutine d_rot2_turnover
