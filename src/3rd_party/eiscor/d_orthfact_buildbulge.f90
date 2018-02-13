#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first Givens' rotations that start a real 
! Francis iteration. If JOB = 'S' a single Givens' rotation is 
! constructed. If JOB = 'D' two Givens' rotations are constructed for
! a real double shift step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DBL             LOGICAL
!                    .TRUE.: real doubleshift
!                    .FALSE.: real singleshift
!
!  Q               REAL(8) array of dimension (4)
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2)
!                    array of generators for complex diagonal matrix
!
!  SHFT            REAL(8) array of dimension (2)
!                    real and imaginary part of complex shift
!
! OUTPUT VARIABLES:
!
!  B1, B2          REAL(8) array of dimension (2)
!                    generators for givens rotations that represent
!                    the first transformations
!                    if DBL == .FALSE. then B2 is unused
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_buildbulge(DBL,Q,D,SHFT,B1,B2)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: DBL
  real(8), intent(in) :: Q(4), D(2), SHFT(2)
  real(8), intent(inout) :: B1(2), B2(2)
  
  ! compute variables
  real(8) :: nrm, COL(3), T(3,2)
  
  
  ! double shift
  if (DBL) then
  
    ! compute first two columns of A
    call d_orthfact_2x2diagblock(.TRUE.,Q,D,T(1:2,1:2))  
    T(3,1) = 0d0
    T(3,2) = Q(4)*D(2)

    ! compute p(A)e_1
    COL(1) = T(1,1)*T(1,1) + T(1,2)*T(2,1) + SHFT(1)*SHFT(1) + SHFT(2)*SHFT(2) - 2d0*T(1,1)*SHFT(1)
    COL(2) = T(2,1)*(T(1,1)+T(2,2)-2d0*SHFT(1))
    COL(3) = T(2,1)*T(3,2)
    
    ! compute first givens rotations
    call d_rot2_vec2gen(COL(2),COL(3),B1(1),B1(2),nrm)
    
    ! compute second givens rotations
    call d_rot2_vec2gen(COL(1),nrm,B2(1),B2(2),nrm)
    
  ! single shift
  else
  
    ! compute first block of A
    call d_orthfact_2x2diagblock(.TRUE.,Q,D,T(1:2,1:2))  

    ! compute first Givens' rotation
    call d_rot2_vec2gen(T(1,1)-SHFT(1),T(2,1),B1(1),B1(2),nrm)
    
  end if

end subroutine d_orthfact_buildbulge

