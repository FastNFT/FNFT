#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_vec4gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generators for a complex rotation represented
! by 3 real numbers: the real and imaginary parts of a complex cosine,
! CR and CI, and a scrictly real sine, S. The CR, CI and S are 
! constructed to be parallel with the vector [AR+iAI,BR+iBI]^T, where AR, 
! AI, BR and BI are real and i = sqrt(-1).
!
! This routine always adjusts the phase of [AR+iAI,BR+iBI]^T so that S >= 0.
!
! If any part of the input vector [AR+iAI,BR+iBI]^T contains a NAN then
! CR, CI, S and NRM are all set to NAN.
!
! If only one of AR, AI, BR or BI = +/-INF then the corresponding AR, AI, BR 
! or BI is first set to +/-1 and the remaining terms are set to 0. Then the 
! CR, CI and S are computed from the new vector containing +/-1 and 0. In 
! this case NRM is always set to INF.
!
! If more than one of AR, AI, BR or BI = +/- INF then CR = CI = S = NRM = NAN
!
! If AR = AI = BR = BI = 0 then CR = 1, CI = S = 0 and NRM = 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EXCEPTIONAL CASES
!
!    AR   |    AI   |    BR   |    BI   |    CR   |    CI   |    S    |   NRM 
! ------- | ------- | ------- | ------- | ------- | ------- | ------- | -------
!   0d0   |   0d0   |   0d0   |   0d0   |   1d0   |   0d0   |   0d0   |   0d0
! ------- | ------- | ------- | ------- | ------- | ------- | ------- | -------
! +-INF   |   XdX   |   XdX   |   XdX   | +-1d0   |   0d0   |   0d0   |   INF
!   XdX   | +-INF   |   XdX   |   XdX   |   0d0   | +-1d0   |   0d0   |   INF
!   XdX   |   XdX   | +-INF   |   XdX   |   0d0   |   0d0   |   1d0   |   INF
!   XdX   |   XdX   |   XdX   | +-INF   |   0d0   |   0d0   |   1d0   |   INF
!          at least two +-INFs          |   NAN   |   NAN   |   NAN   |   NAN
! ------- | ------- | ------- | ------- | ------- | ------- | ------- | -------
!           at least one NAN            |   NAN   |   NAN   |   NAN   |   NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  AR, AI          REAL(8) 
!                    real and imaginary part of the first component 
!                    of the complex vector [AR+iAI,BR+iBI]^T
!
!  BR, BI          REAL(8) 
!                    the second component 
!                    of the complex vector [AR+iAI,BR+iBI]^T
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
!                    of the complex vector [AR+iAI,BR+iBI]^T
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)

  implicit none
  
  ! input variables
  real(8), intent(in) :: AR, AI, BR, BI
  real(8), intent(inout) :: CR, CI, S, NRM
  
  ! compute variables
  real(8), parameter :: inf = EISCOR_DBL_INF
  real(8) :: pAr, pAi, pBr, pBi
  
  ! compute phase of BR, BI
  call d_rot2_vec2gen(BR,BI,pBr,pBi,S)
  
  ! compute phase of AR, AI
  call d_rot2_vec2gen(AR,AI,pAr,pAi,CR)  
  
  ! check if AR is the only INF 
  if ((abs(S)<=inf).AND.(abs(CR)>inf).AND.(pAr.EQ.0d0)) then
  
    ! construct CR, CI, S
    call z_rot3_vec3gen(AR,AI,S,CR,CI,S,NRM)  
      
  ! check if AI is the only INF 
  else if ((abs(S)<=inf).AND.(abs(CR)>inf).AND.(pAi.EQ.0d0)) then
  
    ! construct CR, CI, S
    call z_rot3_vec3gen(AR,AI,S,CR,CI,S,NRM) 
    
  ! otherwise
  else
  
    ! adjust phase of AR, AI and BR, BI so that BR = sqrt(|BR|^2 + |BI|^2) and BI = 0 
    call d_rot2_vec2gen(pAr*pBr + pAi*pBi,-pAr*pBi + pAi*pBr,pAr,pAi,CI)  
   
    ! construct CR, CI, S
    call z_rot3_vec3gen(CR*pAr,CR*pAi,S,CR,CI,S,NRM)
    
  end if

end subroutine z_rot3_vec4gen
