#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_chasedown 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chases the a single misfit core transformation down 
! one row in a factored unitary plus rank one (upr1fpen) matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
!
!  P               LOGICAL array of dimension (2)
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
!  M               INTEGER
!                    leading dimension of V and W
!
!  V,W             COMPLEX(8) array of dimension (M,2)
!                    right schur vectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update V to store right schurvectors 
!
!  MISFIT          REAL(8) array of dimension (3)
!                    array of generators for misfit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_chasedown(VEC,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,MISFIT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M
  logical, intent(inout) :: P(2)
  real(8), intent(inout) :: Q(6), D1(4), C1(6), B1(6)
  real(8), intent(inout) :: D2(4), C2(6), B2(6), MISFIT(3)
  complex(8), intent(inout) :: V(M,2), W(M,2)
  
  ! compute variables
  real(8) :: G1(3), G2(3), G3(3)
  complex(8) :: A(2,2) 
  
  ! set cores for turnover based on P(1)
  ! hess 
  if (.NOT.P(1)) then
 
    G1 = Q(1:3)
    G2 = Q(4:6)
    G3 = MISFIT

  ! invhess
  else

    G1 = MISFIT
    G2 = Q(4:6)
    G3 = Q(1:3)

  end if

  ! compute turnover
  call z_rot3_turnover(G1,G2,G3)

  ! move misfit to appropriate side based on P(2)
  ! hess
  if (.NOT.P(2)) then
 
    ! update Q
    Q(1:3) = G1
    Q(4:6) = G2

    ! update left schurvectors with G3
    if (VEC) then 
    
      A(1,1) = cmplx(G3(1),G3(2),kind=8)
      A(2,1) = cmplx(G3(3),0d0,kind=8)
      A(1,2) = cmplx(-G3(3),0d0,kind=8)
      A(2,2) = cmplx(G3(1),-G3(2),kind=8)

      W = matmul(W,A)
   
    end if

    ! set MISFIT as inverse of G3
    MISFIT(1) = G3(1)
    MISFIT(2) = -G3(2)
    MISFIT(3) = -G3(3)

    ! pass MISFIT through R2
    call z_upr1utri_rot3swap(.TRUE.,D2,C2,B2,MISFIT)

    ! invert MISFIT
    MISFIT(2) = -MISFIT(2)
    MISFIT(3) = -MISFIT(3)
    
    ! update right schurvectors with MISFIT
    if (VEC) then
    
      A(1,1) = cmplx(MISFIT(1),MISFIT(2),kind=8)
      A(2,1) = cmplx(MISFIT(3),0d0,kind=8)
      A(1,2) = cmplx(-MISFIT(3),0d0,kind=8)
      A(2,2) = cmplx(MISFIT(1),-MISFIT(2),kind=8)

      V = matmul(V,A)
   
    end if

    ! pass MISFIT through R1
    call z_upr1utri_rot3swap(.FALSE.,D1,C1,B1,MISFIT)

  ! invhess
  else

    ! update Q
    Q(1:3) = G1
    Q(4:6) = G3

    ! pass G2 through R1
    call z_upr1utri_rot3swap(.TRUE.,D1,C1,B1,G2)

    ! invert G2
    G2(2) = -G2(2)
    G2(3) = -G2(3)

    ! update right schurvectors using G2
    if (VEC) then
     
      A(1,1) = cmplx(G2(1),G2(2),kind=8)
      A(2,1) = cmplx(G2(3),0d0,kind=8)
      A(1,2) = cmplx(-G2(3),0d0,kind=8)
      A(2,2) = cmplx(G2(1),-G2(2),kind=8)

      V = matmul(V,A)

    end if

    ! pass G2 through R2
    call z_upr1utri_rot3swap(.FALSE.,D2,C2,B2,G2)

    ! update left schurvectors using G2
    if (VEC) then
     
      A(1,1) = cmplx(G2(1),G2(2),kind=8)
      A(2,1) = cmplx(G2(3),0d0,kind=8)
      A(1,2) = cmplx(-G2(3),0d0,kind=8)
      A(2,2) = cmplx(G2(1),-G2(2),kind=8)

      W = matmul(W,A)

    end if

    ! copy inverse of G2 to MISFIT
    MISFIT(1) = G2(1)
    MISFIT(2) = -G2(2)
    MISFIT(3) = -G2(3)

  end if

  ! update position flag
  P(1) = P(2)

end subroutine z_upr1fpen_chasedown
