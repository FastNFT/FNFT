#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1utri_rot3swap
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes a rotation through an unitary plus rank one
! upper-triangular matrix (upr1utri) stored as the product of a diagonal 
! matrix and 2 sequences of rotations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DIR             LOGICAL
!                    .TRUE.: pass rotation from left to right
!                    .FALSE.: pass rotation from right to left
!
!  D               REAL(8) array of dimension (4)
!                    array of generators for complex diagonal matrice
!                    in the upper-triangular matrix
!
!  C               REAL(8) array of dimension (6)
!                    first array of generators for upper-triangular part
!
!  B               REAL(8) array of dimension (6)
!                    second array of generators for upper-triangular part
!
!  G               REAL(8) array of dimension (3)
!                    generators for rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1utri_rot3swap(DIR,D,C,B,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: DIR
  real(8), intent(inout) :: D(4), C(6), B(6), G(3)
 
  ! compute variables
  logical :: SYM
  real(8) :: T(3)

  ! initialize SYM
  if ((C(1).EQ.B(1)).AND.(C(2).EQ.-B(2)).AND.(C(3).EQ.-B(3)).AND. &
             (C(4).EQ.B(4)).AND.(C(5).EQ.-B(5)).AND.(C(6).EQ.-B(6))) then
    SYM = .TRUE.
  else
    SYM = .FALSE.
  end if 

  ! L2R
  if (DIR) then
  
    ! through D
    call z_rot3_swapdiag(D,G)
 
    ! through C and B, symmetric case
    if (SYM) then
    
      ! store G
      T = G

      ! through C
      call z_rot3_turnover(C(1:3),C(4:6),T)

      ! update B
      B(1) = C(1)
      B(2) = -C(2)
      B(3) = -C(3)
      B(4) = C(4)
      B(5) = -C(5)
      B(6) = -C(6)

    ! general case
    else

      ! though C
      call z_rot3_turnover(C(1:3),C(4:6),G)
      
      ! through B
      call z_rot3_turnover(B(4:6),B(1:3),G)
    
    end if
  
  ! R2L
  else
  
    ! through C and B, symmetric case
    if (SYM) then
    
      ! store G
      T = G

      ! through B
      call z_rot3_turnover(B(1:3),B(4:6),T)

      ! update B
      C(1) = B(1)
      C(2) = -B(2)
      C(3) = -B(3)
      C(4) = B(4)
      C(5) = -B(5)
      C(6) = -B(6)

    ! general case
    else
  
      ! through B
      call z_rot3_turnover(B(1:3),B(4:6),G)
      
      ! through C
      call z_rot3_turnover(C(4:6),C(1:3),G)
      
    end if
    
    ! through D
    call z_rot3_swapdiag(D,G)
  
  end if

end subroutine z_upr1utri_rot3swap
