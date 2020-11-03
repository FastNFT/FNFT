#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_1Darray_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a one dimensional double array for INFs and NANs.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    positive integer, dimension of A
!
!  A               REAL(8) array of dimension (N)
!                    real array to be checked
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies A contains valid numbers
!                    .FALSE. implies N < 1 or A contains at least one invalid 
!                            number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_1Darray_check(N,A,FLAG)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(in) :: A(N)
  logical, intent(inout) :: FLAG
  
  ! compute variables
  integer :: ii
  
  ! initialize FLAG
  FLAG = .FALSE.
  
  ! check N
  if (N < 1) then
    return
  end if
  
  ! check array
  do ii = 1,N
  
    ! check for INF or NAN
    call d_scalar_check(A(ii),FLAG)
    if (.NOT.FLAG) then
      return
    end if
    
  end do

end subroutine d_1Darray_check
