#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_scalar_argument
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the argument of a complex number A.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  A              COMPLEX(8)
!                    complex number
!
! OUTPUT VARIABLES:
!
!  ARG            REAL(8)
!                    arg of A in the interval [-pi,pi)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_scalar_argument(A,ARG)

  ! input variables
  complex(8), intent(in) :: A
  real(8), intent(inout) :: ARG
  
  ! compute variables
  real(8), parameter :: PI = 3.14159265358979323846264338327950d0
  real(8) :: Ar, Ai
  integer :: INFO
  
  ! initialize INFO
  INFO = 0  
  
  ! check error in debug mode
  if (DEBUG) then
  
    ! check A
    call z_scalar_check(A,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"A is invalid",INFO,-1)
      return
    end if

 end if

 ! set Ar and Ai
 Ar = dble(A)
 Ai = aimag(A)
 
 ! compute arg 1
 if ((abs(Ar).EQ.0d0).AND.(abs(Ai).EQ.0d0)) then
   ARG = 0d0
   
   ! abs(Ar) > abs(Ai)
 else if (abs(Ar) > abs(Ai)) then
   
   ! compute ARG between [-pi/2,pi/2]
   ARG = atan(abs(Ai/Ar))
   
   ! abs(Ar) < abs(Ai)
 else
   
   ! compute ARG between [-pi/2,pi/2]
   ARG = atan(abs(Ar/Ai))
   ARG = PI/2d0 - ARG
   
 end if
 
 ! correct for [0,2pi)
 ! second quadrant
 if ((Ai >= 0).AND.(Ar < 0)) then
   ARG = PI-ARG
   
   ! third quadrant
 else if ((Ai < 0).AND.(Ar < 0)) then
   ARG = PI+ARG
   
   ! fourth quadrant
 else if ((Ai < 0).AND.(Ar > 0)) then
   ARG = 2d0*PI-ARG
 end if
 
 ! shift by pi
 ARG = ARG - PI
    
end subroutine z_scalar_argument
