#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first transformation for a single shift
! iteration on a factored unitary plus rank one (upr1fpen) matrix pencil.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  P               LOGICAL
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (4)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (6)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
!
!  SHFT            COMPLEX(8) 
!                    contains the shift needed for the first transformation
!
! OUTPUT VARIABLES:
!
!  G               REAL(8) array of dimension (3)
!                    on exit contains the generators for the first
!                    transformation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_buildbulge(P,Q,D1,C1,B1,D2,C2,B2,SHFT,G)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: P
  real(8), intent(in) :: Q(3), D1(4), C1(6), B1(6)
  real(8), intent(in) :: D2(4), C2(6), B2(6)
  complex(8), intent(in) :: SHFT
  real(8), intent(inout) :: G(3)
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: R1(2,2), R2(2,2), vec1(2), vec2(2)
  
  ! get R1
  call z_upr1utri_decompress(.FALSE.,2,D1,C1,B1,R1)
     
  ! get R2
  call z_upr1utri_decompress(.FALSE.,2,D2,C2,B2,R2)
     
  ! compute first columns
  ! P == FALSE
  if (.NOT.P) then
   
    ! first column of R1
    vec1(1) = cmplx(Q(1),Q(2),kind=8)*R1(1,1)
    vec1(2) = cmplx(Q(3),0d0,kind=8)*R1(1,1)
    
    ! first column of R2
    vec2(1) = R2(1,1)
    vec2(2) = R2(2,1)
  
  ! P == TRUE
  else
    
    ! Q*e1
    vec2(1) = cmplx(Q(1),-Q(2),kind=8)
    vec2(2) = cmplx(-Q(3),0d0,kind=8)
    
    ! back solve with R1
    vec2(2) = vec2(2)/R1(2,2)
    vec2(1) = (vec2(1) - R1(1,2)*vec2(2))/R1(1,1)
    
    ! multiply by R2
    vec2 = matmul(R2,vec2)

    ! R2^-1 e1
    vec1(1) = cmplx(1d0,0d0,kind=8)
    vec1(2) = cmplx(0d0,0d0,kind=8)
  
  end if
  
  ! insert shift
  vec1 = vec1 - SHFT*vec2

  ! compute eliminator
  call z_rot3_vec4gen(dble(vec1(1)),aimag(vec1(1)),dble(vec1(2)),&
                      aimag(vec1(2)),G(1),G(2),G(3),nrm)
      
end subroutine z_upr1fpen_buildbulge
