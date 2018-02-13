#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_symtrid_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the real Schur factorization of a real 
! symmetric tridiagonal matrix T.
!
! This routine performs a Cayley transformation of the symmetric
! tridiagonal matrix to a unitary matrix (a descending and a 
! ascending sequence of core transformations). The unitary matrix
! is transformed to upper Hessenberg form by core chasing and then 
! passed to z_unifact_qr.
!
! The Cayley (Moebius) transform -(z-i)/(z+i) maps the real line to
! the unit circle, in particular the interval [-1,1] is mapped to
! [-pi/2, pi/2]. 
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
!  SCA             LOGICAL
!                    .TRUE.: scale and shift the matrix to 
!                            have eigenvalues in [-1,1]
!                    .FALSE.: do not scale nor shift the matrix
!                  !! CAUTION: Not scaling and shifting the matrix can 
!                              result in inaccurate eigenvalues !! 
! 
!  N               INTEGER
!                    dimension of matrix
!
!  D               REAL(8) array of dimension (N)  
!                    diagonal entries of T
!                    on exit: D contains the eigenvalues of T
!
!  E               REAL(8) array of dimension (N)
!                    subdiagonal entries of T
!                    on exit: if SCA=.TRUE., E is scaled 
!
!  WORK            REAL(8) array of dimension (5*N)
!                    work space for eigensolver
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z               COMPLEX(8) array of dimension (M,N)
!                    components of schurvectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. and ID = .TRUE. initializes Z to I 
!                    if VEC = .TRUE. and ID = .FALSE. assumes Z initialized
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 2 implies QR failed to converge
!                    INFO = 1 implies factorization failed
!                    INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_symtrid_qr(VEC,ID,SCA,N,D,E,WORK,M,Z,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID, SCA
  integer, intent(in) :: N, M
  real(8), intent(inout) :: WORK(5*N) ! WORK(1:3N-3) = Q ! WORK(3N+1:5N) = QD
  integer, intent(inout) :: ITS(N-1), INFO
  real(8), intent(inout) :: D(N), E(N-1)
  complex(8), intent(inout) :: Z(M,N)
  ! compute variables
  integer :: ii
  real(8) :: scale
   
  ! initialize INFO
  INFO = 0

  if (VEC.AND.ID) then
    Z = cmplx(0d0,0d0,kind=8)
    do ii=1,min(M,N)
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
   end do
  end if

  ! factorize \Phi(T) and reduction to unitary Hessenberg form
  call d_symtrid_factor(VEC,.FALSE.,SCA,N,D,E,WORK(1:(3*(N-3))),WORK((3*N+1):(5*N)),scale,M,Z,INFO)

  ! check info
  if (INFO.EQ.-56) then
     ! matrix is zero, set eigenvalues to zero
     D = 0d0
     INFO = 0
     return
  end if
  if (INFO.NE.0) then 
     ! print error in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"d_symtrid_factor failed",INFO,INFO)
     end if
     INFO = 1
     ! no eigenvalues found, no back transform necessary
     return
  end if

  ! compute eigenvalues
  call z_unifact_qr(VEC,.FALSE.,N,WORK(1:3*N-3),WORK((3*N+1):(5*N-1)),M,Z,ITS,INFO)

  ! check info
  if (INFO.NE.0) then 
     ! print error in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_unifact_qr failed",INFO,INFO)
     end if
     INFO = 2
     ! since some of the eigenvalues have been found, the back transform is performed for all
     ! return
  end if

  ! back transformation
  do ii=1,N
     D(ii) = WORK(3*N+2*ii)/(1d0+WORK(3*N+2*ii-1))
  end do
  
  ! reverse scaling
  if (SCA) then
     do ii=1,N
        D(ii) = D(ii)*scale
     end do
  end if

end subroutine d_symtrid_qr
