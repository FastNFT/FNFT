!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_poly_roots 
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
!  RESIDUALS       COMPLEX(8) array of dimension (N)
!                    residuals of the computed roots
!
!  INFO            INTEGER 
!                    INFO = 1 implies companion QZ algorithm failed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_poly_roots(N,COEFFS,ROOTS,RESIDUALS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: INFO
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)
  real(8), intent(inout) :: RESIDUALS(N)
  
  ! compute variables
  integer :: ii
  real(8) :: scl
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D1(:),C1(:),B1(:)
  real(8), allocatable :: D2(:),C2(:),B2(:)
  real(8) :: normc
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

  ! compute the norm of the monic polynomial
  normc = 0d0
  do ii=1,N
    normc = normc + abs(COEFFS(ii+1)/COEFFS(1))**2
  end do
  normc = sqrt(normc)

  ! our latest analysis shows that using QR (on the equivalent monic polynomial)
  ! is as accurate as using QZ with abs(COEFFS(1)) < 1, but faster
  ! for very large normc the QZ is still advantegous, hence
  ! we choose conservatievly 1e4
  if (normc.LT.1e8) then
    ! use QR
    
    ! fill V 
    sclc = COEFFS(1)
    V(N) = ((-1d0)**(N))*COEFFS(N+1)/sclc
    do ii=1,(N-1)
      V(ii) = -COEFFS(N+1-ii)/sclc
    end do
    
    ! factor companion matrix
    call z_compmat_compress(N,P,V,Q,D1,C1,B1)
    
    ! call z_upr1fpen_qz
    call z_upr1fact_qr(.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,N,V,ITS,INFO)
    
    if (INFO.NE.0) then
      INFO = 1
    end if

    ! extract roots
    call z_upr1utri_decompress(.TRUE.,N,D1,C1,B1,ROOTS)
    
  else
    ! use QZ

    ! fill V and W
    scl = maxval(abs(COEFFS))
    V(N) = ((-1d0)**(N))*COEFFS(N+1)/scl
    do ii=1,(N-1)
      V(ii) = -COEFFS(N+1-ii)/scl
    end do
    W = cmplx(0d0,0d0,kind=8)
    W(N) = COEFFS(1)/scl
    
    ! factor companion matrix
    call z_comppen_compress(N,P,V,W,Q,D1,C1,B1,D2,C2,B2)
    
    ! call z_upr1fpen_qz
    call z_upr1fpen_qz(.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,N,V,W,ITS,INFO)
    
    if (INFO.NE.0) then
      INFO = 1
    end if

    ! extract roots
    call z_upr1utri_decompress(.TRUE.,N,D1,C1,B1,V)
    call z_upr1utri_decompress(.TRUE.,N,D2,C2,B2,W)
    do ii=1,N
      ROOTS(ii) = V(ii)/W(ii)
    end do

    print*, ROOTS
    
  end if
        
  ! compute residuals
  call z_poly_residuals(N,COEFFS,ROOTS,0,RESIDUALS)
    
  ! free memory
  deallocate(P,ITS,Q,D1,C1,B1,D2,C2,B2,V,W)

end subroutine z_poly_roots
