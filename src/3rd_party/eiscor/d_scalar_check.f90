#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_scalar_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a real number to see if it is +-INF or NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NUM             REAL(8) 
!                    real number to be checked
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies valid real(8)
!                    .FALSE. implies NUM is +-INF or NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_scalar_check(NUM,FLAG)

  implicit none
  
  ! input variables
  real(8), intent(in) :: NUM
  logical, intent(inout) :: FLAG
  
  ! initialize INFO
  FLAG = .FALSE.
  
  ! check for +-INF 
  if (abs(NUM)>EISCOR_DBL_INF) then
    return
  end if
  
  ! check for NAN 
  if (NUM.NE.NUM) then
    return
  end if
  
  ! set FLAG to .TRUE.
  FLAG = .TRUE.

end subroutine d_scalar_check
