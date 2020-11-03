#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2array_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the validity of an array of rot2 generators.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    positive integer, number of rot2 generators in A
!
!  A               REAL(8) array of dimension (2*N)
!                    real array of rot2 generators to be checked
!                    generators are stored sequentially: A(1:4) = {C1,S1,C2,S2}
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies A contains valid rot2 generators
!                    .FALSE. implies N < 1 or A contains at least one invalid 
!                            rot2 generator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2array_check(N,A,FLAG)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(in) :: A(2*N)
  logical, intent(inout) :: FLAG
  
  ! compute variables
  integer :: ii
  
  ! initialize FLAG
  FLAG = .FALSE.
  
  ! check N
  if (N < 1) then
    return
  end if
  
  ! check A
  do ii = 1,N
  
    ! check single generator
    call d_rot2_check(A(2*ii-1),A(2*ii),FLAG)
    if (.NOT.FLAG) then
      return
    end if
    
  end do

end subroutine d_rot2array_check
