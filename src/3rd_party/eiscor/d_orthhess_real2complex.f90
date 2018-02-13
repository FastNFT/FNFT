#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthhess_real2complex
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine converts the real Schur factorization of a real
! orthogonal matrix to the complex Schur factorization.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvectors
!                    .FALSE.: no schurvectors
!
!  N               INTEGER
!                    dimension of matrix
!
!  H               REAL(8) array of dimension (N,N)
!                    block diagonal orthogonal matrix
!
!  M               INTEGER
!                    leading dimension of Z
!
!  Z               REAL(8) array of dimension (M,N)
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. contains real Schur vectors 
!
! OUTPUT VARIABLES:
!
!  E               COMPLEX(8) array of dimension (N)
!                    contains eigenvalues
!
!  V               COMPLEX(8) array of dimension (M,N)
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. contains complex schurvectors 
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -2 implies N is invalid
!                    INFO = -3 implies H is invalid
!                    INFO = -4 implies M is invalid
!                    INFO = -5 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthhess_real2complex(VEC,N,H,M,Z,E,V,INFO)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, M
  integer, intent(inout) :: INFO
  real(8), intent(in) :: H(N,N), Z(M,N)
  complex(8), intent(inout) :: E(N), V(M,N)
  
  ! compute variables
  logical :: flg
  integer :: ii, jj, ind
  complex(8) :: block(2,2) 

  ! initialize INFO
  INFO = 0
 
  ! check factorization
    call d_2Darray_check(N,N,H,flg)
  if (.NOT.flg) then
    INFO = -3
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__ &
      ,"H is invalid",INFO,INFO)
    end if
    return
  end if
  
  ! check M
  if (VEC.AND.(M < 1)) then
    INFO = -4
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"M must be at least 1",INFO,INFO)
    end if
    return
  end if
  
  ! check Z
  if (VEC) then
    call d_2Darray_check(M,N,Z,flg)
    if (.NOT.flg) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
      end if
      return
    end if
  end if   
  
  ! initialize storage
  if (VEC) then

    ! initialize V
    do ii=1,M
      do jj=1,N
        V(ii,jj) = cmplx(Z(ii,jj),0d0,kind=8)
      end do
    end do

    ! initialize block
    block(1,1) = cmplx(0d0,1d0,kind=8)
    block(2,1) = cmplx(1d0,0d0,kind=8)
    block(1,2) = cmplx(-1d0,0d0,kind=8)
    block(2,2) = cmplx(0d0,-1d0,kind=8)
    block = block/sqrt(2d0)

  end if
  
  ! loop through eigenvalues
  ind = 1
  do while (ind.LT.N)

    ! single eigenvalue
    if (H(ind+1,ind).EQ.0d0) then

      ! store eigenvalue
      E(ind) = cmplx(H(ind,ind),0d0,kind=8)
      ind = ind + 1

    ! conjugate pair
    else

      ! store eigenvalues
      E(ind) = cmplx(H(ind,ind),H(ind+1,ind),kind=8)
      E(ind+1) = conjg(E(ind))

      ! update eigenvectors
      if (VEC) then
        V(:,ind:(ind+1)) = matmul(V(:,ind:(ind+1)),block)
      end if

      ! check to see that next rotation is identity
      if ((ind.LT.(N-1)).AND.(H(ind+1,ind+2).NE.0d0)) then
        INFO = -3
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__  &
          ,"Not a valid real Schur form.",INFO,INFO)
        end if
        return
      end if

      ! update ind 
      ind = ind + 2

    end if
    
  end do

  ! update last eigenvalue if not conjugate pair
  if (H(N-1,N).EQ.0d0) then

    ! store eigenvalue
    E(N) = cmplx(H(N,N),0d0,kind=8)

  end if
  
end subroutine d_orthhess_real2complex
