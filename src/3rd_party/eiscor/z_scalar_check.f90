#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_scalar_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a complex number to see if it is +-INF or NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NUM             COMPLEX(8) 
!                    real number to be checked
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies valid complex(8)
!                    .FALSE. implies NUM contains +-INF or NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_scalar_check(NUM,FLAG)

  implicit none
  
  ! input variables
  complex(8), intent(in) :: NUM
  logical, intent(inout) :: FLAG
  
  ! initialize INFO
  FLAG = .FALSE.
  
  ! check real part of NUM
  call d_scalar_check(dble(NUM),FLAG)
  if (.NOT.FLAG) then
    return
  end if
  
  ! check imaginary part of NUM
  call d_scalar_check(aimag(NUM),FLAG)

end subroutine z_scalar_check
