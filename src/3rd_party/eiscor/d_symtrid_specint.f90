#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_symtrid_specint
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes upper and lower bounds on the spectrum of a 
! real symmetric tridiagonal matrix T.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NEWT            LOGICAL
!                    .TRUE.: uses Newton's method to refine boundaries
!                    .FALSE.: only uses Gerschgorin bounds
!
!  N               INTEGER
!                    dimension of matrix
!
!  D               REAL(8) array of dimension (N)  
!                    diagonal entries of T
!
!  E               REAL(8) array of dimension (N-1)
!                    subdiagonal entries of T
!
! OUTPUT VARIABLES:
!
!  A,B             REAL(8)
!                    upper and lower bounds on the spectrum of T
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies A or B is NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_symtrid_specint(NEWT,N,D,E,A,B,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: NEWT
  integer, intent(in) :: N
  real(8), intent(in) :: D(N), E(N-1)
  real(8), intent(inout) :: A, B
  integer, intent(inout) :: INFO

  ! compute variables
  integer :: ii
  real(8) :: ta, tb
  real(8) :: tol = EISCOR_DBL_EPS
  
  ! initialize INFO
  INFO = 0

  ! initialize A and B
  A = D(1) - abs(E(1))
  B = D(1) + abs(E(1))

  ! loop through diagonals
  do ii=1,N-2

    ! compute ta and tb
    ta = D(ii) - abs(E(ii)) - abs(E(ii+1))
    tb = D(ii) + abs(E(ii)) + abs(E(ii+1))

    ! check A
    if (ta < A) then
      A = ta
    end if

    ! check B
    if (tb > B) then
      B = tb
    end if

  end do

  ! check final entry
  ta = D(N) - abs(E(N-1))
  tb = D(N) + abs(E(N-1))
  
  ! check A
  if (ta < A) then
    A = ta
  end if

  ! check B
  if (tb > B) then
    B = tb
  end if

  ! return if A and B are 0
  if ((A.EQ.0d0).AND.(B.EQ.0d0)) then
    return
  end if

  ! Newton updates
  if (NEWT) then

    ! maximum of 10 corrections
    do ii=1,10

      ! correction for A
      call d_symtrid_newtonstep(N,D,E,A,ta,INFO)
       
      ! check if INFO = -1
      if (INFO.EQ.-1) then
        ta = tol
        INFO = 0
      end if
      
      ! correction for B
      call d_symtrid_newtonstep(N,D,E,B,tb,INFO)
      
      ! check if INFO = -1
      if (INFO.EQ.-1) then
        tb = -tol
        INFO = 0
      end if
      
      ! update if not converged
      if ((abs(ta) > tol*abs(A)).AND.(abs(tb) > tol*abs(B))) then

        A = A - ta
        B = B - tb
 
      ! exit otherwise
      else
  
        exit

      end if

    end do

  end if

  ! check if A and B are NAN
  if ((A.NE.A).OR.(B.NE.B)) then 
    INFO = -1
  end if

end subroutine d_symtrid_specint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_symtrid_newtonstep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a Newton correction for a real number x of a
! real symmetric tridiagonal matrix T using Sturm sequences.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  D               REAL(8) array of dimension (N)  
!                    diagonal entries of T
!
!  E               REAL(8) array of dimension (N)
!                    subdiagonal entries of T
!
!  X               REAL(8)
!                    input guess 
!
! OUTPUT VARIABLES:
!
!  F               REAL(8)
!                    computed Newton correction at X
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies zero derivative at X
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_symtrid_newtonstep(N,D,E,X,F,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(in) :: D(N), E(N-1)
  real(8), intent(in) :: X
  real(8), intent(inout) :: F
  integer, intent(inout) :: INFO

  ! compute variables
  integer :: ii
  real(8) :: p0, p1, p2
  real(8) :: d0, d1, d2
  
  ! initialize INFO
  INFO = 0

  ! initialize p0, p1 and p2
  p0 = 1d0 
  p1 = D(1) - X
  p2 = (D(2)-X)*p1 - E(1)*E(1)*p0

  ! initialize d0, d1 and d2
  d0 = 0d0 
  d1 = -1d0
  d2 = -p1 + (D(2)-X)*d1 - E(1)*E(1)*d0

  ! loop through diagonals
  do ii=3,N

    ! update d0, d1 and d2
    d0 = -p2 + (D(ii)-X)*d2 - E(ii-1)*E(ii-1)*d1
    d1 = d2
    d2 = d0

    ! update p0, p1 and p2
    p0 = (D(ii)-X)*p2 - E(ii-1)*E(ii-1)*p1
    p1 = p2
    p2 = p0

  end do

  ! compute F
  ! if derivative is zero set Newton step
  ! to zero and throw error
  if (d2.EQ.0d0) then
     F = 0
     INFO = -1
  else
     F = p2/d2
  end if

end subroutine d_symtrid_newtonstep
