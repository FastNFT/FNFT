#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! u_randomseed_initialize
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine initializes the random number generator using the cpu
! clock.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
! OUTPUT VARIABLES:
!
! INFO            INTEGER
!                    INFO = 1 implies array allocation failed
!                    INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine u_randomseed_initialize(INFO)

  implicit none
  
  ! input variables
  integer, intent(inout) :: INFO

  ! compute variables
  integer :: ii, n, clock
  integer, allocatable :: seed(:)
  
  ! initialize INFO
  INFO = 0
  
  ! get size of see        
  call random_seed(size = n)
  
  ! allocate memory for seed
  allocate(seed(n))
  
  ! check allocation
  if (allocated(seed).EQV..FALSE.) then
    INFO = 1
  
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Array allocation failed",INFO,INFO)
    end if 
    
    return
  end if 
  
  ! get current clock time        
  call system_clock(count=clock)
  
  ! store seeds        
  seed = n + clock + 37 * (/ (ii - 1, ii = 1, n) /)

  ! set the generator
  call random_seed(put = seed)
  
  ! free memory        
  deallocate(seed)

end subroutine u_randomseed_initialize
