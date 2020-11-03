#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine diagonalizes a unitary upper hessenberg matrix that is stored as
! a product of N Givens rotations, without computing square roots.
!
! | u1       -v1 |
! | v1  conj(u1) | | u2       -v2 | 
!                  | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                   | v3  conj(u3) | | 0  u4 |                                   
!                                                                     
! The square root free algorithm only requires the storage of the vi^2,
! so the arrays U and VV contain the following:
!
!  U(i) = ui
! VV(i) = vi^2
!
! The input must satisfy the following:
!
!  |U(i)|^2 + VV(i) = 1
!             VV(N) = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for Givens rotations
!                    on output contains eigenvalues
!
!  VV              REAL(8) array of dimension N
!                    array of real generators (squared) for Givens rotations
!
! OUTPUT VARIABLES:
!
!  ITS             INTEGER array of dimension N-1
!                    contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 1 implies no convergence
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies U or VV is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_qr(N,U,VV,ITS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)
  integer, intent(inout) :: INFO, ITS(N-1)
  
  ! compute variables
  integer :: ii, kk 
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  real(8) :: xx
  complex(8) :: nu
  
  ! initialize info
  INFO = 0
  
  ! check N
  if (N < 2) then
    INFO = -1
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be >= 2.",INFO,INFO)
    end if
    return
  end if
  

  ! check U and VV
  do ii = 1,N
    if (abs(abs(U(ii))**2+VV(ii)-1d0) > 10d0*EISCOR_DBL_EPS) then
      INFO = -2
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"U or VV is invalid.",INFO,INFO)
      end if
      return
    end if
  end do
 
  ! initialize storage
  ITS = 0
  
  ! initialize indices
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = 0

  ! iteration loop
  do kk=1,ITMAX

    ! check for completion
    if(STP <= 0)then    
      ! store eigenvalues in U
      do ii = 1,N-1
        U(N+1-ii) = conjg(U(N-ii))*U(N+1-ii)
        xx = dble(U(N+1-ii))**2 + aimag(U(N+1-ii))**2
        U(N+1-ii) = 5d-1*U(N+1-ii)*(3d0-xx)
      end do
      exit
    end if
    
    ! check for deflation
    call z_urffact_deflationcheck(STP-STR+1,U(STR:STP),VV(STR:STP),ZERO)
    
    if (ZERO.GT.0) then
      ITS(STR+ZERO-1) = ITS(STR+ZERO-1) + ITCNT
      ITCNT = 0
    end if
    
    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
      ! update indices
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if greater than 1x1 chase a bulge
    else

      ! check ZERO
      if (ZERO.GT.0) then
        STR = STR+ZERO
      end if

      ! set nu for top deflations
      if (STR > 1) then
        nu = conjg(U(STR-1))
      else
        nu = cmplx(1d0,0d0,kind=8)
      end if

      ! perform singleshift iteration
      call z_urffact_singlestep(STP-STR+2,U(STR:STP+1),VV(STR:STP+1),nu,ITCNT)
     
      ! update indices
      ITCNT = ITCNT + 1
 
    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if
    
  end do

end subroutine z_urffact_qr
