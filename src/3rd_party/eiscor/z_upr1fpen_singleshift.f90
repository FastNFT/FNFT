#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_singleshift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a shift for a factored unitary plus rank one
! (upr1fpen) matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*N)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (3*N)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
!
! OUTPUT VARIABLES:
!
!  SHFT            COMPLEX(8)
!                    shift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_singleshift(P,Q,D1,C1,B1,D2,C2,B2,SHFT)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: P(2)
  real(8), intent(in) :: Q(6), D1(6), C1(9), B1(9)
  real(8), intent(in) :: D2(6), C2(9), B2(9)
  complex(8), intent(inout) :: SHFT
  
  ! compute variables
  complex(8) :: rho, R1(3,3), R2(3,3), H(2,2), K(2,2)

  ! extract R1
  call z_upr1utri_decompress(.FALSE.,3,D1,C1,B1,R1)
  
  ! extract R2
  call z_upr1utri_decompress(.FALSE.,3,D2,C2,B2,R2)

  ! apply first Q
  if (P(2)) then 

    ! update R2
    H(1,1) = cmplx(Q(4),-Q(5),kind=8)
    H(2,1) = cmplx(-Q(6),0d0,kind=8)
    H(1,2) = cmplx(Q(6),0d0,kind=8)
    H(2,2) = cmplx(Q(4),Q(5),kind=8)
    R2(2:3,:) = matmul(H,R2(2:3,:))

  else

    ! update R1
    H(1,1) = cmplx(Q(4),Q(5),kind=8)
    H(2,1) = cmplx(Q(6),0d0,kind=8)
    H(1,2) = cmplx(-Q(6),0d0,kind=8)
    H(2,2) = cmplx(Q(4),-Q(5),kind=8)
    R1(2:3,:) = matmul(H,R1(2:3,:))

   end if

  ! apply second Q
  if (P(1)) then 

    ! update R2
    H(1,1) = cmplx(Q(1),-Q(2),kind=8)
    H(2,1) = cmplx(-Q(3),0d0,kind=8)
    H(1,2) = cmplx(Q(3),0d0,kind=8)
    H(2,2) = cmplx(Q(1),Q(2),kind=8)
    R2(1:2,:) = matmul(H,R2(1:2,:))

  else

    ! update R1
    H(1,1) = cmplx(Q(1),Q(2),kind=8)
    H(2,1) = cmplx(Q(3),0d0,kind=8)
    H(1,2) = cmplx(-Q(3),0d0,kind=8)
    H(2,2) = cmplx(Q(1),-Q(2),kind=8)
    R1(1:2,:) = matmul(H,R1(1:2,:))

  end if

  ! store ratio of bottom right entries
  rho = R1(3,3)/R2(3,3) 
      
  ! compute eigenvalues and eigenvectors
  call z_2x2array_eig(.TRUE.,R1(2:3,2:3),R2(2:3,2:3),H,K)

  ! wilkinson shift
  if(abs(R1(3,3)/R2(3,3)-rho) < abs(R1(2,2)/R2(2,2)-rho))then
     SHFT = R1(3,3)/R2(3,3)
     ! SHFT = R1(3,3)
     ! SHFT2 = R2(3,3)
  else
     SHFT = R1(2,2)/R2(2,2)
     ! SHFT = R1(2,2)
     ! SHFT2 = R2(2,2)
  end if

  ! avoid INFs and NANs
  if ((SHFT.NE.SHFT).OR.(abs(SHFT) > EISCOR_DBL_INF)) then
    SHFT = cmplx(0d0,0d0,kind=8) ! not sure if this is a good idea?
    !SHFT2 = cmplx(0d0,0d0,kind=8) ! not sure if this is a good idea?
  end if

end subroutine z_upr1fpen_singleshift
