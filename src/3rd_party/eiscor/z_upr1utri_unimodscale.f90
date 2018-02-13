#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1utri_unimodscale
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine scales either a row or column of a unitary plus rank
! one upper-triangular matrix (upr1utri) by a complex unimodular scalar.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ROW             LOGICAL
!                    .TRUE.: scale row
!                    .FALSE.: scale column
!
!  D               REAL(8) array of dimension (2)
!                    array of generators for complex diagonal matrix
!                    in the upper-triangular factor
!
!  C               REAL(8) array of dimension (3)
!                    first array of generators for upper-triangular part
!
!  B               REAL(8) array of dimension (3)
!                    second array of generators for upper-triangular part
!
!  SCL             COMPLEX(8) 
!                    scalar, assumed unimodular
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1utri_unimodscale(ROW,D,C,B,SCL)

  implicit none
  
  ! input variables
  logical, intent(in) :: ROW
  real(8), intent(inout) :: D(2), C(3), B(3)
  complex(8), intent(in) :: SCL
 
  ! compute variables
  real(8) :: nrm
  complex(8) :: temp

  ! update D regardless of row or column
  temp = SCL*cmplx(D(1),D(2),kind=8)
  call d_rot2_vec2gen(dble(temp),aimag(temp),D(1),D(2),nrm)
  
  ! update B and C
  if (.NOT.ROW) then

      ! update B
    temp = SCL*cmplx(B(1),B(2),kind=8)
    call z_rot3_vec3gen(dble(temp),aimag(temp),B(3),B(1),B(2),B(3),nrm)
        
    ! update C
    temp = conjg(SCL)*cmplx(C(1),C(2),kind=8)
    call z_rot3_vec3gen(dble(temp),aimag(temp),C(3),C(1),C(2),C(3),nrm)
  
  end if

end subroutine z_upr1utri_unimodscale
