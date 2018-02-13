#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_compmat_compress
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  COEFFS          COMPLEX(8) array of dimension (N)
!                    coefficients for left triangular factor
!
! OUTPUT VARIABLES:
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C,B             REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_compmat_compress(N,P,COEFFS,Q,D,C,B)       

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  complex(8), intent(inout) :: COEFFS(N)
  
  ! compute variables
  integer :: ii
  real(8) :: phr, phi, nrm, beta
  complex(8) :: temp
  
  ! initialize Q
  Q = 0d0
  do ii = 1,(N-1)
    Q(3*ii) = 1d0
  end do

  ! shuffle COEFFS
  do ii = 1,(N-2)

    ! permute if necessary
    if (P(N-ii-1)) then

      temp = (-1d0)**(N-ii-1)*COEFFS(N-ii)
      COEFFS(2:(N-ii)) = COEFFS(1:(N-ii-1))
      COEFFS(N-ii) = temp    

    end if
 
  end do

  ! compress first triangle
  D = 0d0
  do ii = 1,N
    D(2*ii-1) = 1d0
  end do
  B = 0d0
  C = 0d0

  ! compute the phase of last coefficient
  call d_rot2_vec2gen(dble(COEFFS(N)),aimag(COEFFS(N)),phr,phi,beta)
 
  ! store in D
  D(2*N-1) = phr
  D(2*N) = phi

  ! initialize bottom of C
  call d_rot2_vec2gen(beta,-1d0,C(3*N-2),C(3*N),nrm)

  ! initialize bottom of B
  B(3*N-2) = C(3*N)
  B(3*N) = C(3*N-2)

  ! roll up COEFFS into B and C
  temp = cmplx(nrm,0d0,kind=8)
  do ii = 1,(N-1)

    ! compute new C
    call z_rot3_vec4gen(dble(COEFFS(N-ii)),aimag(COEFFS(N-ii)), & 
                        dble(temp),aimag(temp), & 
                        C(3*(N-ii)-2),C(3*(N-ii)-1),C(3*(N-ii)),nrm)

    ! store new B
    B(3*(N-ii)-2) = C(3*(N-ii)-2)
    B(3*(N-ii)-1) = -C(3*(N-ii)-1)
    B(3*(N-ii)) = -C(3*(N-ii))

    ! update last entry of COEFFS using C
    temp = cmplx(C(3*(N-ii)-2),-C(3*(N-ii)-1),kind=8)*COEFFS(N-ii) &
    + C(3*(N-ii))*temp

  end do

end subroutine z_compmat_compress
