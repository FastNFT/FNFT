#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3array_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the validity of an array of rot3 generators.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    positive integer, number of rot3 generators in A
!
!  A               REAL(8) array of dimension (3*N)
!                    real array of rot3 generators to be checked
!                    generators are stored sequentially: 
!                                  A(1:6) = {CR1,CI1,S1,CR2,CI2,S2}
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies A contains valid rot3 generators
!                    .FALSE. implies N < 1 or A contains at least one invalid 
!                            rot3 generator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3array_check(N,A,FLAG)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(in) :: A(3*N)
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
    call z_rot3_check(A(3*ii-2),A(3*ii-1),A(3*ii),FLAG)
    if (.NOT.FLAG) then
      return
    end if
    
  end do

end subroutine z_rot3array_check
