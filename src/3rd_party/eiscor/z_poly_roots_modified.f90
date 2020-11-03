! This file is a modified version of the original EISCOR file z_poly_roots.f90
!
! Changes:
!   - Residuals are not computed
!   - The roots are no longer printed (forgotten printf?)
!   - Always use QR since it is as good as QZ -> https://arxiv.org/pdf/1611.02435.pdf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_poly_roots_modified
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the roots of a polynomial expressed in the monomial
! basis. The vector of COEFFS is in descending order by degree:
!
! p(x) = COEFFS(1)*x^N + COEFFS(2)*x^N-1 + ... + COEFFS(N)*x + COEFFS(N+1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    degree of the polynomial
!
!  COEFFS          COMPLEX(8) array of dimension (N+1)
!                    coefficients of polynomial ordered from highest
!                    degree coefficient to lowest degree
!
! OUTPUT VARIABLES:
!
!  ROOTS           COMPLEX(8) array of dimension (N)
!                    computed roots
!
!  INFO            INTEGER
!                    INFO = 1 implies companion QR algorithm failed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_poly_roots_modified(N,COEFFS,ROOTS,INFO)

  implicit none

  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: INFO
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)

  ! compute variables
  integer :: ii
  real(8) :: scl
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D1(:),C1(:),B1(:)
  real(8), allocatable :: D2(:),C2(:),B2(:)
  complex(8) :: sclc
  complex(8), allocatable :: V(:),W(:)
  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface
  interface
    function l_upr1fact_inversehess(m,flags)
      logical :: l_upr1fact_inversehess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_inversehess
  end interface
  interface
    function l_upr1fact_cmv(m,flags)
      logical :: l_upr1fact_cmv
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_cmv
  end interface
  interface
    function l_upr1fact_random(m,flags)
      logical :: l_upr1fact_random
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_random
  end interface

  ! allocate memory
  allocate(P(N-2),ITS(N-1),Q(3*(N-1)),D1(2*(N+1)),C1(3*N),B1(3*N))
  allocate(V(N),W(N),D2(2*(N+1)),C2(3*N),B2(3*N))

  ! initialize INFO
  INFO = 0

  ! fill P
  P = .FALSE.

  ! always use QR

  ! fill V
  sclc = COEFFS(1)
  V(N) = ((-1d0)**(N))*COEFFS(N+1)/sclc
  do ii=1,(N-1)
     V(ii) = -COEFFS(N+1-ii)/sclc
  end do

  ! factor companion matrix
  call z_compmat_compress(N,P,V,Q,D1,C1,B1)

  ! call z_upr1fpen_qr
  call z_upr1fact_qr(.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,N,V,ITS,INFO)

  if (INFO.NE.0) then
     INFO = 1
  end if

  ! extract roots
  call z_upr1utri_decompress(.TRUE.,N,D1,C1,B1,ROOTS)

  ! free memory
  deallocate(P,ITS,Q,D1,C1,B1,D2,C2,B2,V,W)

end subroutine z_poly_roots_modified
