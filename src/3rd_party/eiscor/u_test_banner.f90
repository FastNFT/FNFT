#include "eiscor.h"
subroutine u_test_banner(FILENAME) 

  implicit none

  ! input variables
  character(*), intent(in) :: FILENAME
  
  ! compute variables
  integer :: length
  
  ! compute length
  length = len(FILENAME)
  length = length-4

  ! print banner
  write(STDERR,'(a,a,a)',advance='no') FILENAME(1:length),REPEAT(' ',36-length),' ... '
  
end subroutine u_test_banner
