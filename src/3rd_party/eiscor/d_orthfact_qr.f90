#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the real Schur factorization of a real 
! orthogonal upper-Hessenberg matrix that is stored as a product of 
! N-1 Givens rotations and a diagonal matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvectors
!                    .FALSE.: no schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize to Z to identity
!                    .FALSE.: assume Z is already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (N)
!                    array of generators for complex diagonal matrix
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z               REAL(8) array of dimension (M,N)
!                    components of schurvectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. and ID = .TRUE. initializes Z to I 
!                    if VEC = .TRUE. and ID = .FALSE. assumes Z initialized
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 1 implies no convergence
!                    INFO = 0 implies successful computation
!                    INFO = -3 implies N is invalid
!                    INFO = -4 implies Q is invalid
!                    INFO = -5 implies D is invalid
!                    INFO = -6 implies M is invalid
!                    INFO = -7 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_qr(VEC,ID,N,Q,D,M,Z,ITS,INFO)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  integer, intent(inout) :: INFO, ITS(N-1)
  real(8), intent(inout) :: Q(2*(N-1)), D(N), Z(M,N)
  
  ! compute variables
  logical :: flg
  integer :: ii, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  real(8) :: block(2,2) 

  ! initialize INFO
  INFO = 0
 
  ! check factorization
  call d_orthfact_factorcheck(N,Q,D,INFO)
  if (INFO.NE.0) then
    INFO = INFO - 2
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N, Q, or D is invalid",INFO,INFO)
    end if
    return
  end if
  
  ! check M
  if (VEC.AND.(M < 1)) then
    INFO = -6
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"M must be at least 1",INFO,INFO)
    end if
    return
  end if
  
  ! check Z
  if (VEC.AND..NOT.ID) then
    call d_2Darray_check(M,N,Z,flg)
    if (.NOT.flg) then
      INFO = -7
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
      end if
      return
    end if
  end if   
  
  ! initialize storage
  ITS = 0
    
  if (VEC.AND.ID) then
    Z = 0d0
    do ii=1,min(M,N)
      Z(ii,ii) = 1d0
    end do
  end if
  
  ! initialize indices
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = 0

  ! loop for bulgechasing
  do kk=1,ITMAX

    ! check for completion
    if(STP <= 0)then
      exit
    end if

    ! check for deflation
    call d_orthfact_deflationcheck(STP-STR+2,Q((2*STR-1):(2*STP)),D(STR:(STP+1)),ZERO)

    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
    
      ! update indices
      ITS(STR+STP-1) = ITCNT
      ITCNT = 0
      STP = STP - 1
      ZERO = 0
      STR = 1
           
    ! if 2x2 block remove and check again
    else if(STP == (STR+ZERO))then
    
      ! call 2x2 deflation
      call d_orthfact_2x2deflation(VEC,Q((2*STP-1):(2*STP)),D(STP:(STP+1)),M,Z(:,STP:(STP+1)))
    
      ! update indices
      ITS(STR+STP-1) = ITCNT
      ITCNT = 0
      STP = STP - 2
      ZERO = 0
      STR = 1
        
    ! if greater than 2x2 chase a bulge
    else

      ! check STR
      if (ZERO.GT.0) then
        STR = STR+ZERO
      end if

      ! perform singleshift iteration
      call d_orthfact_doublestep(VEC,STP-STR+2,Q((2*STR-1):(2*STP)),D(STR:(STP+1)) &
      ,M,Z(:,STR:(STP+1)),ITCNT)
     
      ! update indices
      ITCNT = ITCNT + 1
 
    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if
    
  end do 
  
end subroutine d_orthfact_qr
