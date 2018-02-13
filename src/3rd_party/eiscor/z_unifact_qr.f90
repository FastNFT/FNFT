#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Schur factorization of
! a unitary upper hessenberg matrix that is stored as a product of 
! N-1 Givens rotations and a complex diagonal matrix. 
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
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for Givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z              COMPLEX(8) array of dimension (M,N)
!                   components of schurvectors
!                   if VEC = .FALSE. unused
!                   if VEC = .TRUE. and ID = .TRUE. initializes Z to I 
!                   if VEC = .TRUE. and ID = .FALSE. assumes Z initialized
!
!  ITS            INTEGER array of dimension (N-1)
!                   contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 1 implies no convergence
!                   INFO = 0 implies successful computation
!                   INFO = -3 implies N is invalid
!                   INFO = -4 implies Q is invalid
!                   INFO = -5 implies D is invalid
!                   INFO = -6 implies M is invalid
!                   INFO = -7 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_qr(VEC,ID,N,Q,D,M,Z,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)
  integer, intent(inout) :: INFO, ITS(N-1)
  complex(8), intent(inout) :: Z(M,N)
  
  ! compute variables
  logical :: flg
  integer :: ii, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  
  ! initialize info
  INFO = 0
  
  ! check factorization
  call z_unifact_factorcheck(N,Q,D,INFO)
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
    call z_2Darray_check(M,N,Z,flg)
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
    Z = cmplx(0d0,0d0,kind=8)
    do ii=1,min(M,N)
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
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
    call z_unifact_deflationcheck(STP-STR+2,Q((3*STR-2):(3*STP)) &
    ,D((2*STR-1):(2*STP+2)),ZERO)
    
    if (ZERO.GT.0) then
      ITS(STR+ZERO-1) = ITS(STR+ZERO-1) + ITCNT
      ITCNT = 0
    end if
    
    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
      ! update indices
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if greater than 1x1 chase a bulge
    else

      ! check ZERO
      if (ZERO.GT.0) then
        STR = STR+ZERO
      end if

      ! perform singleshift iteration
      call z_unifact_singlestep(VEC,STP-STR+2,Q((3*STR-2):(3*STP)),D((2*STR-1):(2*STP+2)) &
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

end subroutine z_unifact_qr
