#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_endchase 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine finaalizes one iteration of Francis' singleshift 
! algorithm for a factored unitary plus rank one (upr1fact) matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE. update schurvectors
!                    .FALSE. no schurvectors
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
!                    leading dimesnion of V 
!
!  V               COMPLEX(8) array of dimension (M,N)
!                    right schur vectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update V to store right schurvectors 
!
!  FLAG            LOGICAL                     
!                    position flag for merging the misfit at the bottom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_endchase(VEC,N,P,Q,D,C,B,M,V,G,FLAG)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, FLAG
  integer, intent(in) :: M, N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N), G(3)
  complex(8), intent(inout) :: V(M,N)
  
  ! compute variables
!  integer :: 
  real(8) :: G1(3), G2(3), G3(3)
  complex(8) :: A(2,2)
  
  ! if N == 2 we do a single fusion
  if ( N.EQ.2 ) then

    ! fuse Q and G, G is now a diagonal rotation
    call z_rot3_fusion(.TRUE.,Q(1:3),G)

    ! G scales the rows of the upper-triangular part
    call z_upr1utri_unimodscale(.TRUE.,D(1:2),C(1:3),B(1:3), &
                                cmplx(G(1),G(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D(3:4),C(4:6),B(4:6), &
                                cmplx(G(1),-G(2),kind=8))

    ! return
    return

  end if
  
  ! set cores for turnover based on P(N-2)
  ! hess 
  if (.NOT.P(N-2)) then
 
    G1 = Q(3*(N-2)-2:3*(N-2))
    G2 = Q(3*(N-1)-2:3*(N-1))
    G3 = G

  ! invhess
  else

    G1 = G
    G2 = Q(3*(N-1)-2:3*(N-1))
    G3 = Q(3*(N-2)-2:3*(N-2))

  end if

  ! compute turnover
  call z_rot3_turnover(G1,G2,G3)

  ! bottom fusion based on FLAG
  ! hess
  if (.NOT.FLAG) then
 
    ! update Q
    Q(3*(N-2)-2:3*(N-2)) = G1
    Q(3*(N-1)-2:3*(N-1)) = G2

    ! update V
    if (VEC) then
     
      A(1,1) = cmplx(G3(1),G3(2),kind=8)
      A(2,1) = cmplx(G3(3),0d0,kind=8)
      A(1,2) = cmplx(-G3(3),0d0,kind=8)
      A(2,2) = cmplx(G3(1),-G3(2),kind=8)

      V(:,N-1:N) = matmul(V(:,N-1:N),A)

    end if

    ! pass G3 through upper-triangular part
    call z_upr1utri_rot3swap(.FALSE.,D(2*N-3:2*N), &
                             C(3*N-5:3*N),B(3*N-5:3*N),G3)

    ! fuse G3 with Q
    call z_rot3_fusion(.TRUE.,Q(3*(N-1)-2:3*(N-1)),G3)

    ! scale rows of upper-triangular part
    call z_upr1utri_unimodscale(.TRUE.,D(2*(N-1)-1:2*(N-1)), &
                                C(3*(N-1)-2:3*(N-1)), & 
                                B(3*(N-1)-2:3*(N-1)), &
                                cmplx(G3(1),G3(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D(2*N-1:2*N), &
                                C(3*N-2:3*N), & 
                                B(3*N-2:3*N), &
                                cmplx(G3(1),-G3(2),kind=8))

  ! invhess
  else

    ! update Q
    Q(3*(N-2)-2:3*(N-2)) = G1
    Q(3*(N-1)-2:3*(N-1)) = G3

    ! pass G2 through upper-triangular part
    call z_upr1utri_rot3swap(.TRUE.,D(2*(N-1)-1:2*N), &
                             C(3*(N-1)-2:3*N),B(3*(N-1)-2:3*N),G2)

    ! update V using G2inv
    if (VEC) then
     
      A(1,1) = cmplx(G2(1),-G2(2),kind=8)
      A(2,1) = cmplx(-G2(3),0d0,kind=8)
      A(1,2) = cmplx(G2(3),0d0,kind=8)
      A(2,2) = cmplx(G2(1),G2(2),kind=8)

      V(:,N-1:N) = matmul(V(:,N-1:N),A)

    end if

    ! fuse G2 with Q
    call z_rot3_fusion(.FALSE.,G2,Q(3*(N-1)-2:3*(N-1)))

    ! move G2 to the otherside and update V
    if (VEC) then
     
      V(:,N-1) = V(:,N-1)*cmplx(G2(1),G2(2),kind=8)
      V(:,N) = V(:,N)*cmplx(G2(1),-G2(2),kind=8)

    end if

    ! scale columns of upper-triangular part
    call z_upr1utri_unimodscale(.FALSE.,D(2*(N-1)-1:2*(N-1)), &
                                C(3*(N-1)-2:3*(N-1)), & 
                                B(3*(N-1)-2:3*(N-1)), &
                                cmplx(G2(1),G2(2),kind=8))
    call z_upr1utri_unimodscale(.FALSE.,D(2*N-1:2*N), &
                                C(3*N-2:3*N), & 
                                B(3*N-2:3*N), &
                                cmplx(G2(1),-G2(2),kind=8))

  end if

  ! update position flag
  if (N.GT.3) then
    P(N-3) = P(N-2)
  end if

  ! update P
  P(N-2) = FLAG

end subroutine z_upr1fact_endchase
