#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_factorcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the factorization input into d_orthfact_qr to make sure
! is represents a unitary hessenberg matrix to machine precision. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                    INFO = 0 implies valid factorization
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies Q is invalid
!                    INFO = -3 implies D is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_factorcheck(N,Q,D,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(2*(N-1)), D(N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  logical :: flg
  integer :: ii
  
  ! initialize INFO
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
  
  ! check Q 
  call d_rot2array_check(N-1,Q,flg)
  if (.NOT.flg) then
    INFO = -2
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
    end if
    return
  end if
  
  ! check D 
  do ii=1,N
    if (abs(D(ii)).NE.1d0) then
      INFO = -3
      ! print error in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
      end if
      return
    end if
  end do
  
end subroutine d_orthfact_factorcheck
