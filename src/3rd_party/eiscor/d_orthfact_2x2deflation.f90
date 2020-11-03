#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_2x2deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine deflates a 2x2 block in an orthogonal 
!  upper-hessenberg matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: eigenvectors, assumes Z already initialized
!                    .FALSE.: no eigenvectors
!
!  Q               REAL(8) array of dimension (2)
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) arrays of dimension (2)
!                    array of generators for complex diagonal matrices
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z               REAL(8) array of dimension (M,2)
!                    eigenvectors
!                    if VEC = .TRUE. updates eigenvectors in Z 
!                    if VEC = .FALSE. unused
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_2x2deflation(VEC,Q,D,M,Z)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  real(8), intent(inout) :: Q(2), D(2)
  integer, intent(in) :: M 
  real(8), intent(inout) :: Z(M,2)
  
  ! compute variables
  real(8) :: nrm, row(2), temp(2,2)
  
  ! D scalar
  if (D(1).EQ.D(2)) then
    
    ! update Q
    Q = sign(1d0,D(1))*Q

    ! update D
    D = 1d0

  ! D not scalar
  else

    ! compute first row of Q*D
    row(1) = Q(1)*D(1)
    row(2) = -Q(2)*D(2)

    ! update D
    D(1) = -sign(1d0,row(1))
    D(2) = -D(1)

    ! update D
    Q(1) = 1d0
    Q(2) = 0d0

    ! subtract off +/-1
    row(1) = row(1) - D(1)

    ! construct eliminator
    call d_rot2_vec2gen(row(2),-row(1),temp(1,1),temp(2,1),nrm)
     
    ! update Z
    if (VEC) then
      temp(1,2) = -temp(2,1)
      temp(2,2) = temp(1,1)
      Z = matmul(Z,temp)
    end if

  end if

end subroutine d_orthfact_2x2deflation
