#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_vec3gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generators for a complex rotation represented
! by 3 real numbers: the real and imaginary parts of a complex cosine,
! CR and CI, and a strictly real sine, S. The CR, CI and S are 
! constructed to be parallel with the vector [AR+iAI,B]^T, where AR, 
! AI and B are real and i = sqrt(-1).
!
! If AR = AI = B = 0 then CR = 1, CI = S = 0 and NRM = 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  AR, AI          REAL(8) 
!                    real and imaginary part of the first component 
!                    of the complex vector [AR+iAI,B]^T
!
!  B               REAL(8) 
!                    the second component 
!                    of the complex vector [AR+iAI,B]^T
!
! OUTPUT VARIABLES:
!
!  CR, CI          REAL(8)
!                    on exit CR = AR/NRM, CI = AI/NRM
!
!  S               REAL(8)
!                    on exit S = B/NRM
!
!  NRM             REAL(8)
!                    on exit contains the norm 
!                    of the complex vector [AR+iAI,B]^T
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)

  implicit none
  
  ! input variables
  real(8), intent(in) :: AR, AI, B
  real(8), intent(inout) :: CR, CI, S, NRM
  
  ! compute variables
  real(8) :: tar, tai, tb
  real(8) :: nar, nai, nb

  ! set local variables
  nar = abs(AR)
  nai = abs(AI)
  nb = abs(B)

  ! AR = AI = B = 0
  if(nar.EQ.0d0 .AND. nai.EQ.0d0 .AND. nb.EQ.0d0)then
  
    CR = 1d0
    CI = 0d0
    S = 0d0
    NRM = 0d0
    
  ! |AR| >= |B| and |AR| >= |AI| 
  else if(nar >= nb .AND. nar >= nai)then
  
    tb = B/AR
    tai = AI/AR
    NRM = sign(sqrt(1d0 + tb*tb + tai*tai),AR)
    CR = 1d0/NRM
    CI = tai*CR
    S = tb*CR
    NRM = AR*NRM
    
  ! |AI| >= |B| and |AI| >= |AR| 
  else if(nai >= nb .AND. nai >= nar)then
  
    tb = B/AI
    tar = AR/AI
    NRM = sign(sqrt(1d0 + tb*tb + tar*tar),AI)
    CI = 1d0/NRM
    CR = tar*CI
    S = tb*CI
    NRM = AI*NRM
    
  ! |B| >= |AR| and |B| >= |AI|     
  else
    tar = AR/B
    tai = AI/B
    NRM = sign(sqrt(1d0 + tai*tai + tar*tar),B)
    S = 1d0/NRM
    CR = tar*S
    CI = tai*S
    NRM = B*NRM
    
  end if

end subroutine z_rot3_vec3gen
