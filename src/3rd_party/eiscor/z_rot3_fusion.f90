#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_fusion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine uses the generators for two Givens rotations represented
! by 3 real numbers: the real and imaginary parts of a complex cosine
! and a scrictly real sine and performs a fusion. !
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  FLAG             LOGICAL 
!                    .TRUE.: store product in G1 and diagonal in G2
!                    .FALSE.: store product in G2 and diagonal in G1
!
!  G1, G2           REAL(8) arrays of dimension (3)
!                    generators for givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_fusion(FLAG,G1,G2)

  implicit none

  ! input variables
  logical, intent(in) :: FLAG
  real(8), intent(inout) :: G1(3), G2(3)

  ! compute variables
  real(8) :: c1r, c1i, s1
  real(8) :: c2r, c2i, s2
  real(8) :: c3r, c3i, s3r, s3i
  real(8) :: phr, phi
  real(8) :: nrm
 
  ! retrieve G1  
  c1r = G1(1)
  c1i = G1(2)
  s1 = G1(3)
     
  ! retrieve G2  
  c2r = G2(1)
  c2i = G2(2)
  s2 = G2(3)
    
  ! compute givens product
  c3r = c1r*c2r - c1i*c2i - s1*s2
  c3i = c1r*c2i + c1i*c2r
  s3r = s1*c2r + c1r*s2
  s3i = s1*c2i - c1i*s2
     
  ! compute phase
  call d_rot2_vec2gen(s3r,s3i,phr,phi,nrm)

  ! store product in G1 and diagonal in G2
  if (FLAG) then

      ! update G1
      c2r = c3r*phr + c3i*phi
      c2i = -c3r*phi + c3i*phr
      s2 = s3r*phr + s3i*phi
      call z_rot3_vec3gen(c2r,c2i,s2,G1(1),G1(2),G1(3),nrm)

      ! set G2
      G2(1) = phr
      G2(2) = phi
      G2(3) = 0d0
  
  ! store product in G2 and diagonal in G1
  else

      ! update G2
      c2r = c3r*phr - c3i*phi
      c2i = c3r*phi + c3i*phr
      s2 = s3r*phr + s3i*phi
      call z_rot3_vec3gen(c2r,c2i,s2,G2(1),G2(2),G2(3),nrm)

      ! set G2
      G1(1) = phr
      G1(2) = -phi
      G1(3) = 0d0

  end if
  
end subroutine z_rot3_fusion
