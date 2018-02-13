#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_factorcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the factorization of a factored unitary plus
! rank one (upr1fact) matrix to make sure it represents an extended 
! hessenberg triangular matrix to machine precision. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) arrays of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    in the upper-triangular factor
!
!  C,B             REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular part
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies valid factorization
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies Q is invalid
!                   INFO = -3 implies D is invalid
!                   INFO = -4 implies C is invalid
!                   INFO = -5 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_factorcheck(N,Q,D,C,B,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(in) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  logical :: flg
  integer :: ii
  
  ! initialize INFO
  INFO = 0
  
  ! check N
  if (N < 2) then
    INFO = -1
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if
    return
  end if
  
  ! check Q 
  call z_rot3array_check(N-1,Q,flg)
  if (.NOT.flg) then
    INFO = -2
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
    end if
    return
  end if
  
  ! check D 
  call d_rot2array_check(N,D,flg)
  if (.NOT.flg) then
    INFO = -3
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
    end if
    return
  end if

  ! check C
  call z_rot3array_check(N,C,flg)
  if (.NOT.flg) then
    INFO = -4
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"C is invalid",INFO,INFO)
    end if
    return
  end if

  ! check C for diagonal rotations
  do ii=1,N
    if (C(3*ii).EQ.0d0) then
      INFO = -4
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__, &
             "C contains a diagonal rotation",INFO,INFO)
      end if
      return
   end if
  end do

  ! check B
  call z_rot3array_check(N,B,flg)
  if (.NOT.flg) then
    INFO = -5
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"B is invalid",INFO,INFO)
    end if
    return
  end if

  ! check B for diagonal rotations
  do ii=1,N
    if (B(3*ii).EQ.0d0) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__, &
             "B contains a diagonal rotation",INFO,INFO)
      end if
      return
   end if
  end do

end subroutine z_upr1fact_factorcheck
