#include "eiscor.h"
subroutine u_test_passed(TIME) 

  implicit none

  ! input variables
  real(8), intent(in) :: TIME

  ! print failure
  write(STDERR,'(a,F10.5,a)') 'PASSED in ',TIME,' secs'
  
end subroutine u_test_passed
