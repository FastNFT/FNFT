#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_decompress
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine decompresses a factored unitary plus rank one 
! (upr1fpen) matrix pencil.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
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
!  H,T             COMPLEX(8) arrays of dimension (N,N)
!                    extended hessenberg, triangular pencil
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_decompress(N,P,Q,D1,C1,B1,D2,C2,B2,H,T)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  real(8), intent(in) :: Q(3*(N-1)), D1(2*N), C1(3*N), B1(3*N)
  real(8), intent(in) :: D2(2*N), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: H(N,N), T(N,N)
  
  ! decompress H
  call z_upr1fact_decompress(N,P,Q,D1,C1,B1,H)

  ! decompress T
  call z_upr1utri_decompress(.FALSE.,N,D2,C2,B2,T)

end subroutine z_upr1fpen_decompress
