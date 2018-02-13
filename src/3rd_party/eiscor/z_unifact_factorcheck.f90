#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_factorcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the factorization input into z_unifact_qr to make 
! sure it represents a unitary hessenberg matrix to machine precision. 
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
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies valid factorization
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies Q is invalid
!                   INFO = -3 implies D is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_factorcheck(N,Q,D,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(in) :: Q(3*N-3), D(2*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  logical :: flg
  
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
  call z_rot3array_check(N-1,Q,flg)
  if (.NOT.flg) then
    INFO = -2
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
    end if
    return
  end if
  
  ! check D 
  call d_rot2array_check(N,D,flg)
  if (.NOT.flg) then
    INFO = -3
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
    end if
    return
  end if
  
end subroutine z_unifact_factorcheck
