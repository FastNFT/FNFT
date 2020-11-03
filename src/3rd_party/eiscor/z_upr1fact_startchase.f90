#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_startchase 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine initializes one iteration of Francis' singleshift 
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
! OUTPUT VARIABLES:
!
!  G               REAL(8) array of dimension 3
!                    generators for bulge core transformation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_startchase(VEC,N,P,Q,D,C,B,M,V,ITCNT,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M, N, ITCNT
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N), G(3)
  complex(8), intent(inout) :: V(M,N)
  
  ! compute variables
  integer :: ir1, ir2, id1, id2
  logical :: tp(2)
  real(8) :: Ginv(3)
  real(8) :: tq(6), td(6), tc(9), tb(9)
  complex(8) :: shift, A(2,2)
  
  ! compute shift
  ! random shift
  if ((mod(ITCNT,20).EQ.0).AND.(ITCNT.GT.0)) then
    call random_number(G(1))
    call random_number(G(2))
    shift = cmplx(G(1),G(2),kind=8)
          
  ! wilkinson shift
  else
  
    ! special case N = 2
    if (N.LT.3) then 

      ! pad with identity
      tp = .FALSE.
      tq = 0d0; tq(1) = 1d0; tq(4:6) = Q
      td = 0d0; td(1) = 1d0; td(3:6) = D
      tc = 0d0; tc(3) = 1d0; tc(4:9) = C
      tb = 0d0; tb(3) = -1d0; tb(4:9) = B
    
    ! general case
    else 

      ! store in temp arrays    
      if (N.EQ.3) then
        tp(1) = .FALSE.
        tp(2) = P(N-2)
      else
        tp = P((N-3):(N-2))
      end if
      ir2 = 3*N; ir1 = ir2-8
      id2 = 2*N; id1 = id2-5
      tq = Q((ir1):(ir2-3))
      td = D(id1:id2)
      tc = C(ir1:ir2)
      tb = B(ir1:ir2)

    end if

    ! compute wilkinson shift
    call z_upr1fact_singleshift(tp,tq,td,tc,tb,shift)

  end if

!print*,""
!print*," shift =",shift
!print*,""

  ! build bulge
  call z_upr1fact_buildbulge(P(1),Q(1:3),D(1:4),C(1:6),B(1:6),shift,G)

  ! set Ginv
  Ginv(1) = G(1)
  Ginv(2) = -G(2)
  Ginv(3) = -G(3)
  
  ! update V
  if (VEC) then
    
    A(1,1) = cmplx(G(1),G(2),kind=8)
    A(2,1) = cmplx(G(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    V(:,1:2) = matmul(V(:,1:2),A)
    
  end if
  
  ! initialize turnover 
  ! hess
  if (.NOT.P(1)) then
  
    ! fuse Ginv and Q, Ginv is now a diagonal rotation
    call z_rot3_fusion(.FALSE.,Ginv,Q(1:3))

    ! pass G through triangular part
    call z_upr1utri_rot3swap(.FALSE.,D(1:4),C(1:6),B(1:6),G)
  
    ! move Ginv to the other side and update V
    if (VEC) then
     
      V(:,1) = V(:,1)*cmplx(Ginv(1),Ginv(2),kind=8)
      V(:,2) = V(:,2)*cmplx(Ginv(1),-Ginv(2),kind=8)
     
    end if
  
    ! Ginv scales the columns of the upper-triangular part
    call z_upr1utri_unimodscale(.FALSE.,D(1:2),C(1:3),B(1:3), &
                                cmplx(Ginv(1),Ginv(2),kind=8))
    call z_upr1utri_unimodscale(.FALSE.,D(3:4),C(4:6),B(4:6), &
                                cmplx(Ginv(1),-Ginv(2),kind=8))

  ! inverse hess
  else
  
    ! pass G through triangular part
    call z_upr1utri_rot3swap(.FALSE.,D(1:4),C(1:6),B(1:6),G)
  
    ! fuse Q and G, G is now a diagonal rotation
    call z_rot3_fusion(.TRUE.,Q(1:3),G)

    ! G scales the rows of the upper-triangular part
    call z_upr1utri_unimodscale(.TRUE.,D(1:2),C(1:3),B(1:3), &
                                cmplx(G(1),G(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D(3:4),C(4:6),B(4:6), &
                                cmplx(G(1),-G(2),kind=8))

    ! replace G with Ginv as misfit
    G = Ginv

  end if
  
end subroutine z_upr1fact_startchase
