#include "eiscor.h"
subroutine u_test_failed(LINENUM) 

  implicit none

  ! input variables
  integer, intent(in) :: LINENUM

  ! print failure
  write(STDERR,'(a,I4)') 'FAILED on line: ',LINENUM
  stop
  
end subroutine u_test_failed
