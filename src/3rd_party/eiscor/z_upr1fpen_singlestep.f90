#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a factored unitary plus rank one (upr1fpen) matrix.. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
!
!  FUN             LOGICAL FUNCTION FUN(N,P)
!                    takes integer N and logical array P of 
!                    dimension N-2 and outputs a logical 
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*N)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (3*N)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
!
!  M               INTEGER
!                    leading dimension of V and W
!
!  V,W             COMPLEX(8) array of dimension (M,N)
!                    right and left schurvectors 
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_singlestep(VEC,FUN,N,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M, N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*N), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(M,N),W(M,N)
  integer, intent(in) :: ITCNT
  interface
    function FUN(m,flags)
      logical :: FUN
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function FUN
  end interface
  
  ! compute variables
  integer :: ii, ir1, ir2, id1, id2
  logical :: final_flag
  real(8) :: MISFIT(3)
  
  ! compute final_flag
  if (N.LT.3) then 
    final_flag = .FALSE.
  else
    final_flag = FUN(N,P)
  end if  

  ! initialize core chasing
  call z_upr1fpen_startchase(VEC,N,P,Q,D1,C1,B1,D2,C2,B2,M,V(:,1:2),W(:,1:2),ITCNT,MISFIT)
  
  ! core chasing loop
  do ii=1,(N-3)

    ! compute indices
    ir1 = 3*(ii)+1
    ir2 = 3*(ii+2)
    id1 = 2*(ii)+1
    id2 = 2*(ii+2)

    ! move misfit down one row
    call z_upr1fpen_chasedown(VEC,P(ii:ii+1),Q((ir1-3):(ir2-3)),D1(id1:id2), & 
                              C1(ir1:ir2),B1(ir1:ir2),D2(id1:id2),C2(ir1:ir2), &
                              B2(ir1:ir2),M,V(:,ii+1:ii+2),W(:,ii+1:ii+2),MISFIT)

  end do

  ! finish core chasing
  call z_upr1fpen_endchase(VEC,N,P,Q,D1,C1,B1,D2,C2,B2,M,V(:,N-1:N),W(:,N-1:N),MISFIT,final_flag)
  
end subroutine z_upr1fpen_singlestep
