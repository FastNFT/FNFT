#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthhess_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the real Schur factorization of a real 
! orthogonal upper-Hessenberg matrix H.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvectors
!                    .FALSE.: no schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize to Z to identity
!                    .FALSE.: assume Z is already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  H               REAL(8) array of dimension (N,N)
!                    orthogonal hessenberg matrix, assumed that 
!                    H(ii,jj) = 0 for |ii-jj| > 0
!                    on exit H is a block diagonal matrix
!
!  WORK            REAL(8) array of dimension (3*N)
!                    work space for eigensolver
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z               REAL(8) array of dimension (M,N)
!                    components of schurvectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. and ID = .TRUE. initializes Z to I 
!                    if VEC = .TRUE. and ID = .FALSE. assumes Z initialized
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 2 implies d_orthfact_qr failed
!                    INFO = 1 implies d_orthhess_factor failed
!                    INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthhess_qr(VEC,ID,N,H,WORK,M,Z,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  real(8), intent(inout) :: WORK(3*N)
  integer, intent(inout) :: ITS(N-1), INFO
  real(8), intent(inout) :: H(N,N), Z(M,N)
  
  ! compute variables
  integer :: ii
  
  ! initialize INFO
  INFO = 0
  
  ! compress H
  call d_orthhess_factor(N,H,WORK(1:(2*N)),WORK((2*N+1):(3*N)),INFO)

  ! check info
  if (INFO.NE.0) then 
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_orthhess_factor failed" &
      ,INFO,INFO)
    end if 
    INFO = 1
    return
  end if
  
  ! compute eigenvalues
  call d_orthfact_qr(VEC,ID,N,WORK(1:(2*N)),WORK((2*N+1):(3*N)),M,Z,ITS,INFO) 

  ! check info
  if (INFO.NE.0) then 
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_orthfact_qr failed",INFO,INFO)
    end if
    INFO = 2 
    return
  end if
  
  ! update H
  H = 0d0
  do ii=1,N
    H(ii,ii) = 1d0
  end do
  do ii=1,N-1
    if (abs(WORK(2*ii)).NE.0d0) then 
       H(ii,ii) = WORK(2*ii-1)
       H(ii+1,ii+1) = WORK(2*ii-1)
       H(ii+1,ii) = WORK(2*ii)
       H(ii,ii+1) = -WORK(2*ii)
    end if
  end do
  do ii=1,N
    H(:,ii) = H(:,ii)*WORK(2*N+ii)
  end do
  
end subroutine d_orthhess_qr
