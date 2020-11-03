!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1utri_decompress
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine decompresses a unitary plus rank one upper-triangular 
! matrix (upr1utri). Such a matrix is stored as a product 
! of two sequences of Givens' rotations and a complex unimodular 
! diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DIAG            CHARACTER
!                    .TRUE.: diagonals only
!                    .FALSE.: entire triangular part
!
!  N               INTEGER
!                    dimension of matrix
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  C               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of rotations
!
!  B               REAL(8) array of dimension (3*N)
!                    array of generators for third sequence of rotations
!
! OUTPUT VARIABLES:
!
!  T               COMPLEX(8) array of dimension (N,N)
!                    if DIAG == .TRUE. diagonal of T stored in first column
!                    if DIAG == .FALSE. entire triangular part is constructed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1utri_decompress(DIAG,N,D,C,B,T)

  implicit none
  
  ! input variables
  logical, intent(in) :: DIAG
  integer, intent(in) :: N
  real(8), intent(in) :: D(2*N), C(3*N), B(3*N)
  complex(8), intent(inout) :: T(N,N)
  
  ! compute variables
  integer :: ii,ind
  complex(8) :: g, p, temp(2,2)
  
  ! diagonals only
  if (DIAG) then
   
    ! initialize T
    T(:,1) = cmplx(0d0,0d0,kind=8)

    ! compute diagonals
    do ii = 1,N
      T(ii,1) = -cmplx(D(2*ii-1),D(2*ii),kind=8)*B(3*ii)/C(3*ii)
    end do

  ! entire upper-triangular part
  else

    ! initialize T
    T = cmplx(0d0,0d0,kind=8)
    do ii=1,N
      T(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  
    ! fill with B
    do ii=1,(N-1)
      ind = 3*(ii-1)
      temp(1,1) = cmplx(B(ind+1),B(ind+2),kind=8)
      temp(2,1) = cmplx(B(ind+3),0d0,kind=8)
      temp(1,2) = -temp(2,1)
      temp(2,2) = conjg(temp(1,1))
      T(1:(ii+1),ii:(ii+1)) = matmul(T(1:(ii+1),ii:(ii+1)),temp)
    end do
    T(:,N) = T(:,N)*cmplx(B(3*(N-1)+1),B(3*(N-1)+2),kind=8)
    
    ! update first row
    g = cmplx(1d0,0d0,kind=8)
    p = g
    do ii=1,(N-1)
      ind = 3*(ii-1)
      T(1,ii) = T(1,ii)-(g*cmplx(B(ind+1),B(ind+2),kind=8)+&
                cmplx(C(ind+1),-C(ind+2),kind=8)*B(ind+3)/C(ind+3))/p
      p = p*C(ind+3)
      g = cmplx(B(ind+1),-B(ind+2),kind=8)*cmplx(C(ind+1),-C(ind+2),kind=8)-&
          g*B(ind+3)*C(ind+3)
    end do
    ind = 3*(N-1)
    T(1,N) = T(1,N)-(g*cmplx(B(ind+1),B(ind+2),kind=8)+&
             cmplx(C(ind+1),-C(ind+2),kind=8)*B(ind+3)/C(ind+3))/p
    
    ! apply C
    do ii=1,(N-1)
      ind = 3*(ii-1)
      temp(1,1) = cmplx(C(ind+1),C(ind+2),kind=8)
      temp(2,1) = cmplx(C(ind+3),0d0,kind=8)
      temp(1,2) = -temp(2,1)
      temp(2,2) = conjg(temp(1,1))
      T(ii:(ii+1),(ii+1):N) = matmul(temp,T(ii:(ii+1),(ii+1):N))
      T(ii,ii) = temp(1,1)*T(ii,ii) + temp(1,2)*T(ii+1,ii)
      T(ii+1,ii) = cmplx(0d0,0d0,kind=8)
    end do
    T(N,N) = T(N,N)*cmplx(C(3*(N-1)+1),C(3*(N-1)+2),kind=8) + &
             cmplx(B(3*(N-1)+3),0d0,kind=8)*cmplx(-C(3*(N-1)+3),0d0,kind=8)
    
    ! apply D
    do ii=1,N
      ind = 2*(ii-1)
      T(ii,:) = T(ii,:)*cmplx(D(ind+1),D(ind+2),kind=8)
    end do
  
  end if  
  
end subroutine z_upr1utri_decompress
