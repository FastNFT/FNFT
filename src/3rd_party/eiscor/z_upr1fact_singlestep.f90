#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a factored unitary plus rank one (upr1fact) matrix.. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE. update schurvectors
!                    .FALSE. no schurvectors
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
!                    array of generators for givens rotations
!
!  D               REAL(8) arrays of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  C,B             REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular part
!
!  M               INTEGER
!                    leading dimesnion of V and W
!
!  V               COMPLEX(8) array of dimension (M,N)
!                    right schur vectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update V to store right schurvectors 
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_singlestep(VEC,FUN,N,P,Q,D,C,B,M,V,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M, N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  complex(8), intent(inout) :: V(M,N)
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
  call z_upr1fact_startchase(VEC,N,P,Q,D,C,B,M,V(:,1:2),ITCNT,MISFIT)
  
  ! core chasing loop
  do ii=1,(N-3)

    ! compute indices
    ir1 = 3*(ii)+1
    ir2 = 3*(ii+2)
    id1 = 2*(ii)+1
    id2 = 2*(ii+2)

    ! move misfit down one row
    call z_upr1fact_chasedown(VEC,P(ii:ii+1),Q((ir1-3):(ir2-3)),D(id1:id2), & 
                              C(ir1:ir2),B(ir1:ir2),M,V(:,ii+1:ii+2),MISFIT)

  end do

  ! finish core chasing
  call z_upr1fact_endchase(VEC,N,P,Q,D,C,B,M,V,MISFIT,final_flag)
  
end subroutine z_upr1fact_singlestep
