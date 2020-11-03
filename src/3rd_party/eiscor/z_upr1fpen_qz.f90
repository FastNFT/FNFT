#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_qz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generalized schur decomposition of an extended
! hessenberg triangular matrix pencil where both matrices in the pencil and
! compressed unitary plus rank one matrices (upr1).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize V and W to I
!                    .FALSE.: assume V and W already initialized
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
! OUTPUT VARIABLES:
!
!  V,W             COMPLEX(8) array of dimension (M,N)
!                    right and left schurvectors 
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 2 implies no convergence 
!                    INFO = 1 random seed initialization failed
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N, Q, D1, C1, B1, D2, C2, or B2 is invalid
!                    INFO = -14 implies V is invalid
!                    INFO = -15 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_qz(VEC,ID,FUN,N,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*N), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(M,N), W(M,N)
  integer, intent(inout) :: INFO, ITS(N-1)
  interface
    function FUN(m,flags)
      logical :: FUN
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function FUN
  end interface
  
  ! compute variables
  logical :: flg
  integer :: ii, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  
  ! initialize info
  INFO = 0
  
  ! initialize random seed
  !call u_randomseed_initialize(INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__, & 
           "Failed to initialize random seed",INFO,INFO)
    end if
    INFO = 1
    return
  end if
  
  ! check factorization
  call z_upr1fact_factorcheck(N,Q,D1,C1,B1,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__, & 
           "N, Q, D1, C1 or B1 is invalid",INFO,INFO)
    end if
    INFO = -1
    return
  end if

  call z_upr1fact_factorcheck(N,Q,D2,C2,B2,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__, & 
           "N, Q, D2, C2 or B2 is invalid",INFO,INFO)
    end if
    INFO = -1
    return
  end if

  ! check V
  if (VEC.AND..NOT.ID) then
    call z_2Darray_check(M,N,V,flg)
    if (.NOT.flg) then
      INFO = -14
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"V is invalid",INFO,INFO)
      end if
      return
    end if
  end if

  ! check W
  if (VEC.AND..NOT.ID) then
    call z_2Darray_check(M,N,W,flg)
    if (.NOT.flg) then
      INFO = -15
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"W is invalid",INFO,INFO)
      end if
      return
    end if
  end if

  ! initialize storage
  ITS = 0
  
  if (VEC.AND.ID) then
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,min(m,n)
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
  
  if (VEC.AND.ID) then
    W = cmplx(0d0,0d0,kind=8)
    do ii=1,min(m,n)
      W(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
  
  ! initialize indices
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = 0
  
  ! iteration loop
  do kk=1,ITMAX

    ! check for completion
    if(STP <= 0)then    
      exit
    end if

    ! check for deflation
    call z_upr1fpen_deflationcheck(VEC,STP-STR+2,P(STR:(STP-1)), &
         Q((3*STR-2):(3*STP)),D1((2*STR-1):(2*STP+2)),C1((3*STR-2):(3*STP+3)), & 
         B1((3*STR-2):(3*STP+3)),D2((2*STR-1):(2*STP+2)),C2((3*STR-2):(3*STP+3)), & 
         B2((3*STR-2):(3*STP+3)),M,V(:,STR:STP+1),W(:,STR:STP+1),ZERO)
    
    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
    
      ! update indices
      ITS(STR+ZERO-1) = ITCNT
      ITCNT = 0
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if greater than 1x1 chase a bulge
    else

      ! update STR
      if (ZERO.GT.0) then
         STR = STR + ZERO
         ZERO = 0
         ITS(STR+ZERO-1) = ITCNT
         ITCNT = 0
      end if

      ! perform singleshift iteration
      call z_upr1fpen_singlestep(VEC,FUN,STP-STR+2,P(STR:(STP-1)),Q((3*STR-2):(3*STP)) &
      ,D1((2*STR-1):(2*STP+2)),C1((3*STR-2):(3*STP+3)),B1((3*STR-2):(3*STP+3)) &
      ,D2((2*STR-1):(2*STP+2)),C2((3*STR-2):(3*STP+3)),B2((3*STR-2):(3*STP+3)) &
      ,M,V(:,STR:(STP+1)),W(:,STR:(STR+1)),ITCNT)
     
      ! update indices
      if (ITCNT.EQ.-1) then 
        ITCNT = 1 
      else
        ITCNT = ITCNT + 1
      end if

    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if
    
  end do

end subroutine z_upr1fpen_qz
