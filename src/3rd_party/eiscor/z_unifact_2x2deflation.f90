#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_2x2deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine deflates a 2x2 block in a unitary upper-hessenberg matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: eigenvectors, assumes Z already initialized
!                    .FALSE.: no eigenvectors
!
!  Q               REAL(8) array of dimension (3)
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) arrays of dimension (4)
!                    array of generators for complex diagonal matrices
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z               COMPLEX(8) array of dimension (M,2)
!                    eigenvectors
!                    if VEC = .TRUE. updates eigenvectors in Z 
!                    if VEC = .FALSE. unused
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_2x2deflation(VEC,Q,D,M,Z)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  real(8), intent(inout) :: Q(3), D(4)
  integer, intent(in) :: M 
  complex(8), intent(inout) :: Z(M,2)
  
  ! compute variables
  real(8) :: nrm, G1(3), G2(3), G3(3), Qt(6)
  complex(8) :: A(2,2), Zt(2,2)
  
  ! compute 2x2 blocks
  Qt = 0d0; Qt(1:3) = Q; Qt(4) = 1d0
  call z_unifact_2x2diagblock(.TRUE.,Qt,D,A)

  !  schur decomposition
  call z_2x2array_eig(.FALSE.,A,A,Zt,Zt)

  ! set Q
  Q(1) = 1d0
  Q(2) = 0d0
  Q(3) = 0d0

  ! set D
  call d_rot2_vec2gen(dble(A(1,1)),aimag(A(1,1)),D(1),D(2),nrm)
  call d_rot2_vec2gen(dble(A(2,2)),aimag(A(2,2)),D(3),D(4),nrm)
    
  ! update Z
  if (VEC) then
    Z = matmul(Z,Zt)
  end if

end subroutine z_unifact_2x2deflation
