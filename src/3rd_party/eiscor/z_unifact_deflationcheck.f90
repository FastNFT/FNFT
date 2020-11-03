#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in a unitary upper hessenberg 
! matrix that is stored as a product of givens rotations and a complex 
! diagonal matrix. When a deflation occurs the corresponding rotation
! is set to the identity matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    generators must be orthogonal to working precision
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     index of the last known deflation
!                     on output contains index of newest deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_deflationcheck(N,Q,D,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ZERO
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)

  ! compute variables
  integer :: ii, jj, down
  real(8), parameter :: tol = EISCOR_DBL_EPS
  real(8) :: dr, di, qr, qi, s, nrm

  ! intialize ZERO
  ZERO = 0
  
  ! check for deflation
  do ii=1,(N-1)
  
    ! deflate if subdiagonal is small enough
    nrm = abs(Q(3*(N-ii)))
    if(nrm < tol)then
        
      ! set ZERO
      ZERO = max(0,N-ii)  ! why 0?

      ! extract diagonal
      qr = Q(3*ZERO-2)
      qi = Q(3*ZERO-1)

      ! set rotation to identity
      Q(3*ZERO-2) = 1d0
      Q(3*ZERO-1) = 0d0
      Q(3*ZERO) = 0d0
        
      ! update first diagonal
      dr = D(2*ZERO-1)
      di = D(2*ZERO)
        
      nrm = qr*dr - qi*di
      di = qr*di + qi*dr
      dr = nrm
      call d_rot2_vec2gen(dr,di,D(2*ZERO-1),D(2*ZERO),nrm)
        
      ! deflate downward
      do jj = 1,(N-1-ZERO)
     
        ! set downward index
        down = ZERO+jj

        ! update Q
        dr = Q(3*down-2)
        di = Q(3*down-1)
        s = Q(3*down)
              
        nrm = qr*dr + qi*di
        di = qr*di - qi*dr
        dr = nrm

        call z_rot3_vec3gen(dr,di,s,Q(3*down-2),Q(3*down-1),Q(3*down),nrm) 
      
      end do
           
      ! update second diagonal
      dr = D(2*N-1)
      di = D(2*N)
           
      nrm = qr*dr + qi*di
      di = qr*di - qi*dr
      dr = nrm
      call d_rot2_vec2gen(dr,di,D(2*N-1),D(2*N),nrm)
        
      ! exit loop  
      exit

    end if

  end do
  
end subroutine z_unifact_deflationcheck
