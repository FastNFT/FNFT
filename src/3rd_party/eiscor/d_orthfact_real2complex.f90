#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_real2complex
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine converts the real Schur factorization of a real 
! orthogonal matrix to the complex Schur factorization. The real 
! orthogonal matrix is stored as the product of N-1 Givens rotations 
! and a diagonal matrix with entries +/- 1.
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
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (N)
!                    array of generators for complex diagonal matrix
!
!  M               INTEGER
!                    leading dimension of Z
!
!  Z               REAL(8) array of dimension (M,N)
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. contains real schurvectors 
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
!                    INFO = -3 implies Q is invalid
!                    INFO = -4 implies D is invalid
!                    INFO = -5 implies M is invalid
!                    INFO = -6 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_real2complex(VEC,N,Q,D,M,Z,E,V,INFO)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, M
  integer, intent(inout) :: INFO
  real(8), intent(in) :: Q(2*(N-1)), D(N), Z(M,N)
  complex(8), intent(inout) :: E(N), V(M,N)
  
  ! compute variables
  logical :: flg
  integer :: ii, jj, ind
  complex(8) :: block(2,2) 

  ! initialize INFO
  INFO = 0
 
  ! check factorization
  call d_orthfact_factorcheck(N,Q,D,INFO)
  if (INFO.NE.0) then
    INFO = INFO - 1
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__ &
      ,"N, Q, or D is invalid",INFO,INFO)
    end if
    return
  end if
  
  ! check M
  if (VEC.AND.(M < 1)) then
    INFO = -5
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
      INFO = -6
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
    if (Q(2*ind).EQ.0d0) then

      ! store eigenvalue
      E(ind) = cmplx(D(ind),0d0,kind=8)
      ind = ind + 1

    ! conjugate pair
    else

      ! store eigenvalues
      E(ind) = cmplx(Q(2*ind-1),Q(2*ind),kind=8)
      E(ind+1) = conjg(E(ind))

      ! update eigenvectors
      if (VEC) then
        V(:,ind:(ind+1)) = matmul(V(:,ind:(ind+1)),block)
      end if

      ! check to see that next rotation is identity
      if (ind.LT.(N-1)) then
         if (Q(2*ind+2).NE.0d0) then
            INFO = -3
            ! print error message in debug mode
            if (DEBUG) then
               call u_infocode_check(__FILE__,__LINE__  &
                    ,"Not a valid real Schur form.",INFO,INFO)
            end if
            return
         end if
      end if

      ! update ind 
      ind = ind + 2

    end if
    
  end do

  ! update last eigenvalue if not conjugate pair
  if (Q(2*(N-1)).EQ.0d0) then

    ! store eigenvalue
    E(N) = cmplx(D(N),0d0,kind=8)

  end if
  
end subroutine d_orthfact_real2complex
