#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine uses the generators for three Givens rotations represented
! by 3 real numbers: the real and imaginary parts of a complex cosine
! and a scrictly real sine and performs a turnover. The input arrays should be 
! ordered:
!
!    G1  G3
!      G2
!
! The new generators are stored as:
!
!      G1
!    G3  G2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  G1, G2, G3       REAL(8) arrays of dimension (3)
!                    generators for givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_turnover(G1,G2,G3)

  implicit none

  ! input variables
  real(8), intent(inout) :: G1(3), G2(3), G3(3)

  ! compute variables
  real(8) :: nrm, ar, ai, b
  real(8) :: c1r, c1i, s1
  real(8) :: c2r, c2i, s2
  real(8) :: c3r, c3i, s3
  real(8) :: c4r, c4i, s4
  real(8) :: c5r, c5i, s5
  real(8) :: c6r, c6i, s6
  
  ! set local variables
  c1r = G1(1)
  c1i = G1(2)
  s1 = G1(3)
  c2r = G2(1)
  c2i = G2(2)
  s2 = G2(3)
  c3r = G3(1)
  c3i = G3(2)
  s3 = G3(3)

  ! compute first rotation
  ar = s1*c3r + (c1r*c2r + c1i*c2i)*s3
  ai = s1*c3i + (c1r*c2i - c1i*c2r)*s3
  b = s2*s3
  call z_rot3_vec3gen(ar,ai,b,c4r,c4i,s4,nrm)
 
  ! nrm == 0
  if (nrm.EQ.0d0) then 
    ! compute second rotation
    ar = c1r*c3r - c1i*c3i - s1*c2r*s3
    ai = c1r*c3i + c1i*c3r - s1*c2i*s3
    call d_rot2_unitvec2gen(ar,ai,c5r,c5i)
    s5 = 0d0

    ! update first rotation
    call d_rot2_vec2gen(c1r,c1i,ar,ai,nrm)
    c4r = c5r*ar + c5i*ai
    c4i = c5i*ar - c5r*ai
    call d_rot2_unitvec2gen(c4r,c4i,c4r,c4i)
    
    ! compute third rotation
    ar = c2r*c4r + c2i*c4i
    ai = c2i*c4r - c2r*c4i
    b = nrm*s2
    call z_rot3_unitvec3gen(ar,ai,b,c6r,c6i,s6)

  ! all other cases
  else 
    ! compute second rotation
    ar = c1r*c3r - c1i*c3i - s1*c2r*s3
    ai = c1r*c3i + c1i*c3r - s1*c2i*s3
    b = nrm
    call z_rot3_unitvec3gen(ar,ai,b,c5r,c5i,s5)
    
    ! compute third rotation
    ar = c2r*c4r + c2i*c4i + c1r*s2*s4
    ai = c2i*c4r - c2r*c4i + c1i*s2*s4
    b = s1*s2/s5
    call z_rot3_unitvec3gen(ar,ai,b,c6r,c6i,s6)
  end if
  
  ! store output
  G1(1) = c5r
  G1(2) = c5i
  G1(3) = s5
  G2(1) = c6r
  G2(2) = c6i
  G2(3) = s6
  G3(1) = c4r
  G3(2) = c4i
  G3(3) = s4
  
end subroutine z_rot3_turnover
