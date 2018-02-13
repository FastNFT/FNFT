#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_startchase 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine initializes one iteration of Francis' singleshift 
! algorithm for a factored unitary plus rank one (upr1fpen) matrix pencil.
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
!  V,W             COMPLEX(8) array of dimension (M,N)
!                    right and left schurvectors 
!
! OUTPUT VARIABLES:
!
!  G               REAL(8) array of dimension 3
!                    generators for bulge core transformation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_startchase(VEC,N,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ITCNT,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M, N, ITCNT
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*N), C2(3*N), B2(3*N), G(3)
  complex(8), intent(inout) :: V(M,2), W(M,2)
  
  ! compute variables
  integer ::  ir1, ir2, id1, id2
  logical :: tp(2)
  real(8) :: Ginv(3)
  real(8) :: tq(6), td1(6), tc1(9), tb1(9)
  real(8) :: td2(6), tc2(9), tb2(9)
  complex(8) :: shift, A(2,2)
  
  ! compute shift
  ! random shift
  if ((mod(ITCNT,15).EQ.0).AND.(ITCNT.GT.0)) then
    call random_number(G(1))
    call random_number(G(2))
    shift = cmplx(G(1),G(2),kind=8)
          
  ! wilkinson shift
  else
  
    ! special case N = 2
    if (N.LT.3) then 

      ! pad with identity
      tp = .FALSE.
      tq = 0d0; tq(1) =  1d0; tq(4:6) = Q
      td1 = 0d0; td1(1) =  1d0; td1(3:6) = D1
      tc1 = 0d0; tc1(3) =  1d0; tc1(4:9) = C1
      tb1 = 0d0; tb1(3) = -1d0; tb1(4:9) = B1
      td2 = 0d0; td2(1) =  1d0; td2(3:6) = D2
      tc2 = 0d0; tc2(3) =  1d0; tc2(4:9) = C2
      tb2 = 0d0; tb2(3) = -1d0; tb2(4:9) = B2
    
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
      td1 = D1(id1:id2)
      tc1 = C1(ir1:ir2)
      tb1 = B1(ir1:ir2)
      td2 = D2(id1:id2)
      tc2 = C2(ir1:ir2)
      tb2 = B2(ir1:ir2)

    end if

    ! compute wilkinson shift
    call z_upr1fpen_singleshift(tp,tq,td1,tc1,tb1,td2,tc2,tb2,shift)

  end if

  ! build bulge
  if (N.LE.2) then
     call z_upr1fpen_buildbulge(.TRUE.,Q(1:3),D1(1:4),C1(1:6),B1(1:6),D2(1:4),C2(1:6),B2(1:6),shift,G)
  else
     call z_upr1fpen_buildbulge(P(1),Q(1:3),D1(1:4),C1(1:6),B1(1:6),D2(1:4),C2(1:6),B2(1:6),shift,G)
  end if

  ! set Ginv
  Ginv(1) = G(1)
  Ginv(2) = -G(2)
  Ginv(3) = -G(3)
  
  ! update left schurvectors with G
  if (VEC) then
    
    A(1,1) = cmplx(G(1),G(2),kind=8)
    A(2,1) = cmplx(G(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    W = matmul(W,A)
    
  end if

  ! copy Ginv into G
  G = Ginv

  ! initialize turnover 
  if (N.LE.2) then

    ! fuse Ginv and Q, Ginv is now a diagonal rotation
    call z_rot3_fusion(.FALSE.,Ginv,Q(1:3))

    ! pass G through R2
    call z_upr1utri_rot3swap(.TRUE.,D2(1:4),C2(1:6),B2(1:6),G)
  
    ! update left schurvectors diagonal Ginv
    if (VEC) then
      
      W(:,1) = W(:,1)*cmplx(Ginv(1),Ginv(2),kind=8)
      W(:,2) = W(:,2)*cmplx(Ginv(1),-Ginv(2),kind=8)
      
    end if

    ! Ginv scales the rows of R2
    call z_upr1utri_unimodscale(.TRUE.,D2(1:2),C2(1:3),B2(1:3), &
                                cmplx(Ginv(1),-Ginv(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D2(3:4),C2(4:6),B2(4:6), &
                                cmplx(Ginv(1),Ginv(2),kind=8))

    ! invert G
    G(2) = -G(2)
    G(3) = -G(3)

    ! update right schurvectors with G
    if (VEC) then
      
      A(1,1) = cmplx(G(1),G(2),kind=8)
      A(2,1) = cmplx(G(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      V = matmul(V,A)
      
    end if

    ! pass G through R1
    call z_upr1utri_rot3swap(.FALSE.,D1(1:4),C1(1:6),B1(1:6),G)

  ! hess
  elseif (.NOT.P(1)) then
  
    ! fuse Ginv and Q, Ginv is now a diagonal rotation
    call z_rot3_fusion(.FALSE.,Ginv,Q(1:3))

    ! pass G through R2
    call z_upr1utri_rot3swap(.TRUE.,D2(1:4),C2(1:6),B2(1:6),G)
  
    ! update left schurvectors diagonal Ginv
    if (VEC) then
      
      W(:,1) = W(:,1)*cmplx(Ginv(1),Ginv(2),kind=8)
      W(:,2) = W(:,2)*cmplx(Ginv(1),-Ginv(2),kind=8)
      
    end if

    ! Ginv scales the rows of R2
    call z_upr1utri_unimodscale(.TRUE.,D2(1:2),C2(1:3),B2(1:3), &
                                cmplx(Ginv(1),-Ginv(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D2(3:4),C2(4:6),B2(4:6), &
                                cmplx(Ginv(1),Ginv(2),kind=8))

    ! invert G
    G(2) = -G(2)
    G(3) = -G(3)

    ! update right schurvectors with G
    if (VEC) then
      
      A(1,1) = cmplx(G(1),G(2),kind=8)
      A(2,1) = cmplx(G(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      V = matmul(V,A)
      
    end if

    ! pass G through R1
    call z_upr1utri_rot3swap(.FALSE.,D1(1:4),C1(1:6),B1(1:6),G)

  ! inverse hess
  else
  
    ! pass Ginv through R2
    call z_upr1utri_rot3swap(.TRUE.,D2(1:4),C2(1:6),B2(1:6),Ginv)
  
    ! invert Ginv
    Ginv(2) = -Ginv(2)
    Ginv(3) = -Ginv(3)

    ! update right schurvectors with Ginv
    if (VEC) then
      
      A(1,1) = cmplx(Ginv(1),Ginv(2),kind=8)
      A(2,1) = cmplx(Ginv(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      V = matmul(V,A)
      
    end if

    ! pass Ginv through R1
    call z_upr1utri_rot3swap(.FALSE.,D1(1:4),C1(1:6),B1(1:6),Ginv)

    ! fuse Ginv and Q, Ginv is now a diagonal rotation
    call z_rot3_fusion(.TRUE.,Q(1:3),Ginv)

    ! Ginv scales the rows of R1
    call z_upr1utri_unimodscale(.TRUE.,D1(1:2),C1(1:3),B1(1:3), &
                                cmplx(Ginv(1),Ginv(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D1(3:4),C1(4:6),B1(4:6), &
                                cmplx(Ginv(1),-Ginv(2),kind=8))

  end if
  
end subroutine z_upr1fpen_startchase
