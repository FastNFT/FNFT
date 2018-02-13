#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the generators for a complex rotation represented
! by 3 real numbers: the real and imaginary parts of a complex cosine,
! CR and CI, and a scrictly real sine, S. 
!
! To check that it is valid, a new rotation XR, XI, Y is computed. If
! |CR-XR|, |CI-XI|, |S-Y|, |NRM-1| < 3*EISCOR_DBL_EPS the generators are 
! considered valid. Otherwise they are invalid.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  CR, CI          REAL(8)
!                    generators for complex cosine CR + iCI
!
!  S               REAL(8)
!                    generator for real sine S
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies valid generators
!                    .FALSE. implies invalid generators
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_check(CR,CI,S,FLAG)

  implicit none
  
  ! input variables
  real(8), intent(in) :: CR, CI, S
  logical, intent(inout) :: FLAG
  
  ! compute variables
  real(8), parameter :: tol = 3d0*EISCOR_DBL_EPS
  real(8) :: XR, XI, Y, NRM

  ! initialize FLAG to .FALSE.
  FLAG = .FALSE.
  
  ! compute new rotation from input one
  call z_rot3_vec3gen(CR,CI,S,XR,XI,Y,NRM)
  
  ! check for equality
  if ((abs(CR-XR)<tol).AND.(abs(CI-XI)<tol).AND.(abs(S-Y)<tol).AND.(abs(1d0-NRM)<tol)) then
    FLAG = .TRUE.
  end if

end subroutine z_rot3_check
