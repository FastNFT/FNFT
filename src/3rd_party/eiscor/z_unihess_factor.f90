#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unihess_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine factors a unitary upper Hessenberg matrix into a 
! sequence of Givens' rotations and a diagonal matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  H               COMPLEX(8) array of dimension (N,N)
!                    unitary matrix to be reduced
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for Givens' rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies H is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unihess_factor(N,H,Q,D,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)
  integer, intent(inout) :: INFO
  complex(8), intent(inout) :: H(N,N)
  
  ! compute variables
  logical :: flg
  integer :: ii
  real(8) :: nrm, tol 
  real(8) :: cr, ci, s 
  
  ! initialize info
  INFO = 0
  
  ! check N
  if (N < 2) then
    INFO = -1
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if
    return
  end if
  
  ! check H
  call z_2Darray_check(N,N,H,flg)
  if (.NOT.flg) then
    INFO = -2
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"H is invalid",INFO,INFO)
    end if
    return
  end if  
  
  ! set tol
  tol = max(10d0,dble(N))*EISCOR_DBL_EPS
  
  ! loop for reduction
  ! ii is the column being reduced
  do ii=1,(N-1)
            
    ! reduce to block diagonal
    call z_rot3_vec4gen(dble(H(ii,ii)),aimag(H(ii,ii)),dble(H(ii+1,ii)),aimag(H(ii+1,ii)),cr,ci,s,nrm)
    
    ! check for unitarity       
    if (abs(nrm-1d0) >= tol) then
      INFO = -2
      ! print error in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"H is not unitary",INFO,INFO)
      end if 
      return          
    end if
    
    ! update H
    H(ii,ii) = cmplx(cr,-ci,kind=8)*H(ii,ii) + cmplx(s,0d0,kind=8)*H(ii+1,ii)
    H(ii+1,(ii+1):N) = cmplx(-s,0d0,kind=8)*H(ii,(ii+1):N) + cmplx(cr,ci,kind=8)*H(ii+1,(ii+1):N)
          
    ! store in Q
    Q(3*(ii-1)+1) = cr
    Q(3*(ii-1)+2) = ci
    Q(3*(ii-1)+3) = s
          
    ! store in D
    call d_rot2_vec2gen(dble(H(ii,ii)),aimag(H(ii,ii)),D(2*(ii-1)+1),D(2*(ii-1)+2),nrm)
   
  end do
  
  ! store in D
  call d_rot2_vec2gen(dble(H(N,N)),aimag(H(N,N)),D(2*(N-1)+1),D(2*(N-1)+2),nrm)
  
  ! check for unitarity       
  if (abs(nrm-1d0) >= tol) then
    INFO = -2
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"H is not unitary",INFO,INFO)
    end if 
    return          
  end if
          
end subroutine z_unihess_factor
