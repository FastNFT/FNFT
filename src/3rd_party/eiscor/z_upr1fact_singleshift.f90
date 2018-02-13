#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_singleshift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a shift for a factored unitary plus rank one
! (upr1fact) matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  P               LOGICAL array of dimension (2)
!                    position flags for Q
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) arrays of dimension (6)
!                    array of generators for complex diagonal matrix
!                    in the upper-triangular factor
!
!  C,B             REAL(8) arrays of dimension (9)
!                    array of generators for upper-triangular part
!
! OUTPUT VARIABLES:
!
!  SHFT            COMPLEX(8)
!                    shift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_singleshift(P,Q,D,C,B,SHFT)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: P(2)
  real(8), intent(in) :: Q(6), D(6), C(9), B(9)
  complex(8), intent(inout) :: SHFT
  
  ! compute variables
  complex(8) :: rho, R1(3,3), R2(3,3), H(2,2), K(2,2)

  ! extract first triangle
  call z_upr1utri_decompress(.FALSE.,3,D,C,B,R1)
  
  ! extract second triangle
  R2 = cmplx(0d0,0d0,kind=8)
  R2(1,1) = cmplx(1d0,0d0,kind=8)
  R2(2,2) = cmplx(1d0,0d0,kind=8)
  R2(3,3) = cmplx(1d0,0d0,kind=8)

!print*,""
!print*,"R1"
!print*,R1(1,:)
!print*,R1(2,:)
!print*,R1(3,:)
!print*,""
!print*,""
!print*,"R2"
!print*,R2(1,:)
!print*,R2(2,:)
!print*,R2(3,:)
!print*,""

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
  else
    SHFT = R1(2,2)/R2(2,2)
  end if

  ! avoid INFs and NANs
  if ((SHFT.NE.SHFT).OR.(abs(SHFT) > EISCOR_DBL_INF)) then
    SHFT = cmplx(1d9,0d0,kind=8) ! not sure if this is a good idea?
  end if

end subroutine z_upr1fact_singleshift
