#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the generators for a real rotation represented
! by 2 real numbers: a real cosine, C, and a real sine, S. 
!
! To check that it is valid, a new rotation X, Y is computed. If
! |C-X|, |S-Y|, |NRM-1| < 2*EISCOR_DBL_EPS the generators are 
! considered valid. Otherwise they are invalid.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  C               REAL(8)
!                    generators for cosine
!
!  S               REAL(8)
!                    generator for sine
!
! OUTPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies valid generators
!                    .FALSE. implies invalid generators
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_check(C,S,FLAG)

  implicit none
  
  ! input variables
  real(8), intent(in) :: C, S
  logical, intent(inout) :: FLAG
  
  ! compute variables
  real(8), parameter :: tol = 2d0*EISCOR_DBL_EPS
  real(8) :: X, Y, NRM

  ! initialize FLAG to .FALSE.
  FLAG = .FALSE.
  
  ! compute new rotation from input one
  call d_rot2_vec2gen(C,S,X,Y,NRM)
  
  ! check for equality
  if ((abs(C-X)<tol).AND.(abs(S-Y)<tol).AND.(abs(1d0-NRM)<tol)) then
    FLAG = .TRUE.
  end if

end subroutine d_rot2_check
