#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_2Darray_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a two dimensional complex array for INFs and NANs.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  M, N            INTEGER
!                    positive integers, dimensions of A
!
!  A               COMPLEX(8) array of dimension (M,N)
!                    complex array to be checked
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies A contains valid numbers
!                    .FALSE. implies M < 1, N < 1 or A contains at least one 
!                            invalid number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_2Darray_check(M,N,A,FLAG)

  implicit none
  
  ! input variables
  integer, intent(in) :: M, N
  complex(8), intent(in) :: A(M,N)
  logical, intent(inout) :: FLAG
  
  ! compute variables
  integer :: ii, jj
  
  ! initialize FLAG
  FLAG = .FALSE.
  
  ! check M and N
  if ((N < 1).OR.(M < 1)) then
    return
  end if
  
  ! check array
  do ii = 1,M
    do jj = 1,N
  
      ! check for NAN
      call z_scalar_check(A(ii,jj),FLAG)
      if (.NOT.FLAG) then
        return
      end if
    
    end do
  end do

end subroutine z_2Darray_check
