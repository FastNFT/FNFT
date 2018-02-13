#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_1Darray_check 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a one dimensional complex array for INFs and NANs.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    positive integer, dimension of A
!
!  A               COMPLEX(8) array of dimension (N)
!                    complex array to be checked
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies A contains valid numbers
!                    .FALSE. implies N < 1 or A contains at least one invalid 
!                            number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_1Darray_check(N,A,FLAG)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(in) :: A(N)
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
    call z_scalar_check(A(ii),FLAG)
    if (.NOT.FLAG) then
      return
    end if
    
  end do

end subroutine z_1Darray_check
