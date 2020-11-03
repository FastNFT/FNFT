#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in a factored unitary plus rank
! one (upr1fact) matrix. When a deflation occurs the corresponding 
! rotation in the unitary part is set to the identity matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: update eigenvectors
!                    .FALSE.: ignore V matrix
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  C               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of rotations
!
!  B               REAL(8) array of dimension (3*N)
!                    array of generators for second sequence of rotations
!
!  M               INTEGER
!                    leading dimension of V
!
!  V               COMPLEX(8) array of dimension (M,N)
!                    array of eigenvectors
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     on output contains index of newest deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_deflationcheck(VEC,N,P,Q,D,C,B,M,V,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N, M
  logical, intent(in) :: VEC, P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  complex(8), intent(inout) :: V(M,N)
  integer, intent(inout) :: ZERO

  ! compute variables
  integer :: ii, down
  real(8), parameter :: tol = EISCOR_DBL_EPS
  real(8) :: qr, qi, nrm
  
  ! check for deflation
  do ii=1,(N-1)
  
    ! deflate if subdiagonal is small enough
    nrm = abs(Q(3*(N-ii)))
    if(nrm < tol)then

      ! set ZERO
      ZERO = max(0,N-ii) ! why 0?
      
      ! extract diagonal
      qr = Q(3*ZERO-2)
      qi = Q(3*ZERO-1)
                
      ! set rotation to identity
      Q(3*ZERO-2) = 1d0
      Q(3*ZERO-1) = 0d0
      Q(3*ZERO) = 0d0
      
      ! upper left entry of Q(ZERO)
      ! deflation at top
      if ( ZERO.EQ.1 ) then

        ! scale row of upper triangular part
        call z_upr1utri_unimodscale(.TRUE.,D(2*ZERO-1:2*ZERO), &
             C(3*ZERO-2:3*ZERO),B(3*ZERO-2:3*ZERO),cmplx(qr,qi,kind=8))

      ! deflation in the middle
      ! P(ZERO-1) == 0
      elseif ( .NOT.P(ZERO-1) ) then

        ! scale row of upper triangular part
        call z_upr1utri_unimodscale(.TRUE.,D(2*ZERO-1:2*ZERO), &
             C(3*ZERO-2:3*ZERO),B(3*ZERO-2:3*ZERO),cmplx(qr,qi,kind=8))

      ! P(ZERO-1) == 1
      else

        ! scale column of upper triangular part
        call z_upr1utri_unimodscale(.FALSE.,D(2*ZERO-1:2*ZERO), &
             C(3*ZERO-2:3*ZERO),B(3*ZERO-2:3*ZERO),cmplx(qr,qi,kind=8))

        ! update eigenvectors
        if (VEC) then
          V(:,ZERO) = V(:,ZERO)*cmplx(qr,qi,kind=8)
        end if

      end if
      
      ! lower entry of Q(ZERO)
      ! deflation at bottom
      if ( ZERO.EQ.(N-1) ) then

        ! scale row of upper triangular part
        call z_upr1utri_unimodscale(.TRUE.,D(2*ZERO+1:2*ZERO+2), &
             C(3*ZERO+1:3*ZERO+3),B(3*ZERO+1:3*ZERO+3),cmplx(qr,-qi,kind=8))

      ! deflation in the middle
      ! P(ZERO) == 0
      elseif ( .NOT.P(ZERO) ) then

        ! scale column of upper triangular part
        call z_upr1utri_unimodscale(.FALSE.,D(2*ZERO+1:2*ZERO+2), &
             C(3*ZERO+1:3*ZERO+3),B(3*ZERO+1:3*ZERO+3),cmplx(qr,-qi,kind=8))

        ! update eigenvectors
        if (VEC) then
          V(:,ZERO+1) = V(:,ZERO+1)*cmplx(qr,-qi,kind=8)
        end if

      ! P(ZERO) == 1
      else

        ! scale row of upper triangular part
        call z_upr1utri_unimodscale(.TRUE.,D(2*ZERO+1:2*ZERO+2), &
             C(3*ZERO+1:3*ZERO+3),B(3*ZERO+1:3*ZERO+3),cmplx(qr,-qi,kind=8))

      end if

      ! exit loop 
      exit
        
    end if
    
  end do
  
end subroutine z_upr1fact_deflationcheck
