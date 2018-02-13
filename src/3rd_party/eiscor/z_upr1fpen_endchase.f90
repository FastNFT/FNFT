#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_endchase 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine finaalizes one iteration of Francis' singleshift 
! algorithm for a factored unitary plus rank one (upr1fpen) matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
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
!  V,W             COMPLEX(8) array of dimension (M,2)
!                    right and left schurvectors 
!
!  FLAG            LOGICAL                     
!                    position flag for merging the misfit at the bottom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_endchase(VEC,N,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,G,FLAG)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, FLAG
  integer, intent(in) :: M, N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*N), C2(3*N), B2(3*N), G(3)
  complex(8), intent(inout) :: V(M,2),W(M,2)
  
  ! compute variables
  real(8) :: G1(3), G2(3), G3(3)
  complex(8) :: A(2,2) 
  
  ! if N == 2 we do a single fusion
  if ( N.EQ.2 ) then

    ! fuse Q and G, G is now a diagonal rotation
    call z_rot3_fusion(.TRUE.,Q(1:3),G)

    ! G scales the rows of the upper-triangular part
    call z_upr1utri_unimodscale(.TRUE.,D1(1:2),C1(1:3),B1(1:3), &
                                cmplx(G(1),G(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D1(3:4),C1(4:6),B1(4:6), &
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

    ! update left schurvectors with G3
    if (VEC) then 
    
      A(1,1) = cmplx(G3(1),G3(2),kind=8)
      A(2,1) = cmplx(G3(3),0d0,kind=8)
      A(1,2) = cmplx(-G3(3),0d0,kind=8)
      A(2,2) = cmplx(G3(1),-G3(2),kind=8)

      W = matmul(W,A)
   
    end if

    ! invert G3
    G3(2) = -G3(2)
    G3(3) = -G3(3)
 
    ! pass G3 through R2
    call z_upr1utri_rot3swap(.TRUE.,D2(2*N-3:2*N), &
                             C2(3*N-5:3*N),B2(3*N-5:3*N),G3)

    ! invert G3
    G3(2) = -G3(2)
    G3(3) = -G3(3)
 
    ! update right schurvectors with G3
    if (VEC) then
    
      A(1,1) = cmplx(G3(1),G3(2),kind=8)
      A(2,1) = cmplx(G3(3),0d0,kind=8)
      A(1,2) = cmplx(-G3(3),0d0,kind=8)
      A(2,2) = cmplx(G3(1),-G3(2),kind=8)

      V = matmul(V,A)
   
    end if

    ! pass G3 through R1
    call z_upr1utri_rot3swap(.FALSE.,D1(2*N-3:2*N), &
                             C1(3*N-5:3*N),B1(3*N-5:3*N),G3)

    ! fuse G3 with Q
    call z_rot3_fusion(.TRUE.,Q(3*(N-1)-2:3*(N-1)),G3)

    ! scale rows of R1
    call z_upr1utri_unimodscale(.TRUE.,D1(2*(N-1)-1:2*(N-1)), &
                                C1(3*(N-1)-2:3*(N-1)), & 
                                B1(3*(N-1)-2:3*(N-1)), &
                                cmplx(G3(1),G3(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D1(2*N-1:2*N), &
                                C1(3*N-2:3*N), & 
                                B1(3*N-2:3*N), &
                                cmplx(G3(1),-G3(2),kind=8))

  ! invhess
  else

    ! update Q
    Q(3*(N-2)-2:3*(N-2)) = G1
    Q(3*(N-1)-2:3*(N-1)) = G3

    ! pass G2 through R1
    call z_upr1utri_rot3swap(.TRUE.,D1(2*(N-1)-1:2*N), &
                             C1(3*(N-1)-2:3*N),B1(3*(N-1)-2:3*N),G2)

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
    call z_upr1utri_rot3swap(.FALSE.,D2(2*(N-1)-1:2*N), &
                             C2(3*(N-1)-2:3*N),B2(3*(N-1)-2:3*N),G2)

    ! update left schurvectors using G2
    if (VEC) then
     
      A(1,1) = cmplx(G2(1),G2(2),kind=8)
      A(2,1) = cmplx(G2(3),0d0,kind=8)
      A(1,2) = cmplx(-G2(3),0d0,kind=8)
      A(2,2) = cmplx(G2(1),-G2(2),kind=8)

      W = matmul(W,A)

    end if

    ! invert G2
    G2(2) = -G2(2)
    G2(3) = -G2(3)

    ! fuse G2 with Q, G2 is now diagonal
    call z_rot3_fusion(.FALSE.,G2,Q(3*(N-1)-2:3*(N-1)))

    ! update left schurvectors with G2
    if (VEC) then
     
      W(:,1) = W(:,1)*cmplx(G2(1),G2(2),kind=8)
      W(:,2) = W(:,2)*cmplx(G2(1),-G2(2),kind=8)

    end if

    ! scale rows of R2
    call z_upr1utri_unimodscale(.TRUE.,D2(2*(N-1)-1:2*(N-1)), &
                                C2(3*(N-1)-2:3*(N-1)), & 
                                B2(3*(N-1)-2:3*(N-1)), &
                                cmplx(G2(1),-G2(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D2(2*N-1:2*N), &
                                C2(3*N-2:3*N), & 
                                B2(3*N-2:3*N), &
                                cmplx(G2(1),G2(2),kind=8))

  end if

  ! update position flag
  if (N.GT.3) then
    P(N-3) = P(N-2)
  end if

  ! update P
  P(N-2) = FLAG

end subroutine z_upr1fpen_endchase
