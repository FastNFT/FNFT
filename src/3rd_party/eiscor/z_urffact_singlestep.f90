#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift algorithm on a
! unitary upper hessenberg matrix that is stored as a product of givens
! rotations. 
!                                                                               
! | nu  0 | | u1       -v1 |
! |  0  1 | | v1  conj(u1) | | u2       -v2 | 
!                            | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                             | v3  conj(u3) | | 0  u4 |                                   
!                                                                               
! The square root free algorithm only requires the storage of the vi^2,
! so the arrays U and VV contain the following:
!
!  U(i) = ui
! VV(i) = vi^2
!
! The input must satisfy the following:
!
!  |U(i)|^2 + VV(i) = 1
!             VV(N) = 0
!              |NU| = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER 
!                    dimension of matrix, must be >= 2
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for Givens rotations
!
!  VV              REAL(8) array of dimension N
!                    array of real generators for Givens rotations
!
!  NU              COMPLEX(8)
!                    unimodular phase 
!
!  ITCNT           INTEGER array of dimension N-1
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_singlestep(N,U,VV,NU,ITCNT)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ITCNT
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)
  complex(8), intent(in) :: NU
  
  ! compute variables
  integer :: ii
  real(8) :: nn, zz, cc, ss, xx
  complex(8) :: z, w, rho
  complex(8) :: ut
  real(8) :: vvt
  complex(8) :: block(2,2), t1(2,2), t2(2,2)

  ! get 2x2 block
  block(1,1) =  U(N-1)
  block(2,2) =  conjg(U(N-1))
  block(1,2) = -sqrt(VV(N-1))
  block(2,1) =  sqrt(VV(N-1))
  block(:,2) =  block(:,2)*U(N)
  if (N > 2) then
    xx = abs(U(N-2))
    if (xx.GT.0) then
      block(1,:) = conjg(U(N-2))*block(1,:)/xx
    end if
  end if
    
  ! compute eigenvalues and eigenvectors
  t1 = block
  call z_2x2array_eig(.FALSE.,t1,t1,t2,t2)
    
  ! choose wikinson shift
  ! complex abs does not matter here
  if(abs(block(2,2)-t1(1,1)) < abs(block(2,2)-t1(2,2)))then
    rho = t1(1,1)
  else
    rho = t1(2,2)
  end if

  ! compute a nonzero shift
  ! random shift
  xx = abs(rho)
  if (xx == 0) then
    call random_number(xx)
    rho = cmplx(cos(xx),sin(xx),kind=8)
  ! wilkinson shift
  else
    rho = rho/xx
  end if

  ! initialize
  w = -rho
  cc = 1d0
  ss = 0d0

  ! main chasing loop
  do ii=0,(N-1)

    ! set ut and vvt
    ut = NU*U(ii+1)
    vvt = VV(ii+1)

    ! turnover
    call z_rfr3_turnover(w,cc,ss,ut,vvt,rho)

    ! store ut and vvt
    if ( ii > 0 ) then
      U(ii) = conjg(NU)*ut
      VV(ii) = vvt
    end if
    
  end do

end subroutine z_urffact_singlestep
