#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unihess_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Schur factorization of a unitary 
! upper-Hessenberg matrix.
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
!  H               COMPLEX(8) array of dimension (N,N)
!                    unitary hessenberg matrix, assumed that 
!                    H(ii,jj) = 0 for |ii-jj| > 0
!                    on exit contains a diagonal matrix whose entries 
!                    are the eigenvalues of H
!
!  WORK            REAL(8) array of dimension (5*N)
!                    work space for eigensolver
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z              COMPLEX(8) array of dimension (M,N)
!                   components of schurvectors
!                   if VEC = .FALSE. unused
!                   if VEC = .TRUE. and ID = .TRUE. initializes Z to I 
!                   if VEC = .TRUE. and ID = .FALSE. assumes Z initialized
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 2 implies z_unifact_qr failed
!                   INFO = 1 implies z_unihess_factor failed
!                   INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unihess_qr(VEC,ID,N,H,WORK,M,Z,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  real(8), intent(inout) :: WORK(5*N)
  integer, intent(inout) :: ITS(N-1), INFO
  complex(8), intent(inout) :: H(N,N), Z(M,N)
  
  ! compute variables
  integer :: ii
  
  ! initialize INFO
  INFO = 0
  
  ! compress H
  call z_unihess_factor(N,H,WORK(1:(3*N)),WORK((3*N+1):(5*N)),INFO)

  ! check info
  if (INFO.NE.0) then 
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_unihess_factor failed",INFO,INFO)
    end if 
    INFO = 1
    return
  end if
  
  ! compute eigenvalues
  call z_unifact_qr(VEC,ID,N,WORK(1:(3*N)),WORK((3*N+1):(5*N)),M,Z,ITS,INFO) 

  ! check info
  if (INFO.NE.0) then 
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_unifact_qr failed",INFO,INFO)
    end if
    INFO = 2 
    return
  end if
  
  ! update H
  H = cmplx(0d0,0d0,kind=8)
  do ii=1,N
    H(ii,ii) = cmplx(WORK((3*N)+2*(ii-1)+1),WORK((3*N)+2*(ii-1)+2),kind=8)
  end do
  
end subroutine z_unihess_qr
