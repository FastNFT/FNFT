#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_2x2array_eig
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generalized Schur decomposition of a 
! 2x2 matrix pencil (A,B). On exit the pencil (A,B) will be in upper-
! triangular form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  FLAG           LOGICAL
!                    .TRUE. makes no assumptions on B (General)
!                    .FALSE. assumes B = I
!
!  A, B           COMPLEX(8) array of dimension (2,2)
!                   The 2x2 matrix pencil. Upper-triangular on exit.
!                   If FLAG = .FALSE., B is not used.
!
! OUTPUT VARIABLES:
!
!  Q, Z           COMPLEX(8) array of dimension (2,2)
!                   On exit unitary transformations such that Q*(A,B)Z is 
!                   upper-triangular. If FLAG = .FALSE., Z is unused.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_2x2array_eig(FLAG,A,B,Q,Z)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: FLAG
  complex(8), intent(inout) :: A(2,2), B(2,2), Q(2,2), Z(2,2)
  
  ! compute variables
  real(8) :: nrm, cr, ci, s
  complex(8) :: temp(2,2), p0, p1, p2, lambda 
  complex(8) :: trace, detm, disc
  
  ! B not the identity
  if (FLAG) then
  
    ! make B upper-triangular
    ! create eliminator
    call z_rot3_vec4gen(dble(B(1,1)),aimag(B(1,1)),dble(B(2,1)) &
    ,aimag(B(2,1)),cr,ci,s,nrm)
      
    ! store in Q
    Q(1,1) = cmplx(cr,ci,kind=8)
    Q(2,1) = cmplx(s,0d0,kind=8)
    Q(2,2) = conjg(Q(1,1))
    Q(1,2) = -Q(2,1)
      
    ! update B
    B = matmul(conjg(transpose(Q)),B)
    B(2,1) = cmplx(0d0,0d0,kind=8)
      
    ! update A
    A = matmul(conjg(transpose(Q)),A)
      
    ! set Z
    Z(1,1) = cmplx(1d0,0d0,kind=8)
    Z(2,1) = cmplx(0d0,0d0,kind=8)
    Z(2,2) = cmplx(1d0,0d0,kind=8)
    Z(1,2) = cmplx(0d0,0d0,kind=8)
    
    ! check for infinite eigenvalues
    if (abs(B(1,1)).EQ.0d0) then
  
      ! create eliminator
      call z_rot3_vec4gen(dble(A(1,1)),aimag(A(1,1)),dble(A(2,1)) &
      ,aimag(A(2,1)),cr,ci,s,nrm)
      
      ! store in temp
      temp(1,1) = cmplx(cr,ci,kind=8)
      temp(2,1) = cmplx(s,0d0,kind=8)
      temp(2,2) = conjg(temp(1,1))
      temp(1,2) = -temp(2,1)
      
      ! update A
      A = matmul(conjg(transpose(temp)),A)
      A(2,1) = cmplx(0d0,0d0,kind=8)
      
      ! update B
      B = matmul(conjg(transpose(Q)),B)
      
      ! update Q
      Q = matmul(Q,temp)
      
      ! return
      return
  
    ! check for infinite eigenvalues
    else if (abs(B(2,2)).EQ.0d0) then
  
      ! create eliminator
      call z_rot3_vec4gen(dble(A(2,2)),aimag(A(2,2)),dble(A(2,1)) &
      ,aimag(A(2,1)),cr,ci,s,nrm)
      
      ! store in temp
      temp(1,1) = cmplx(cr,-ci,kind=8)
      temp(2,1) = cmplx(s,0d0,kind=8)
      temp(2,2) = conjg(temp(1,1))
      temp(1,2) = -temp(2,1)
      
      ! update A
      A = matmul(A,temp)
      A(2,1) = cmplx(0d0,0d0,kind=8)
      
      ! update B
      B = matmul(B,temp)
      
      ! update Z
      Z = matmul(Z,temp)
      
      ! return
      return
  
    ! compute eigenvalues otherwise
    else 
  
      ! compute polynomial coefficients
      p2 = B(1,1)*B(2,2)
      p1 = A(2,1)*B(1,2) - (A(1,1)*B(2,2)+A(2,2)*B(1,1))
      p0 = A(1,1)*A(2,2) - A(2,1)*A(1,2)  
  
      ! compute intermediate values
      disc = sqrt(p1**2 - 4d0*p2*p0)
  
      ! compute most accurate lambda
      if(abs(-p1+disc) > abs(-p1-disc))then
        lambda = (-p1+disc)/2d0/p2
      else
        lambda = (-p1-disc)/2d0/p2
      end if
  
      ! compute right eigenvector 
      temp = A - lambda*B

      ! create eliminator
      call z_rot3_vec4gen(dble(temp(1,1)),-aimag(temp(1,1)) &
      ,dble(temp(1,2)),-aimag(temp(1,2)),cr,ci,s,nrm)
      
      ! store in temp
      temp(1,1) = cmplx(-s,0d0,kind=8)
      temp(2,1) = cmplx(cr,-ci,kind=8)
      temp(1,2) = cmplx(cr,ci,kind=8)
      temp(2,2) = cmplx(s,0d0,kind=8)
  
      ! update A
      A = matmul(A,temp)
  
      ! update B
      B = matmul(B,temp)
  
      ! update Z
      Z = matmul(Z,temp)
  
      ! return to upper-triangular form  
      ! create eliminator
      call z_rot3_vec4gen(dble(A(1,1)),aimag(A(1,1)),dble(A(2,1)) &
      ,aimag(A(2,1)),cr,ci,s,nrm)
      
      ! store in temp
      temp(1,1) = cmplx(cr,ci,kind=8)
      temp(2,1) = cmplx(s,0d0,kind=8)
      temp(1,2) = cmplx(-s,0d0,kind=8)
      temp(2,2) = cmplx(cr,-ci,kind=8)
  
      ! update A
      A = matmul(conjg(transpose(temp)),A)
      A(2,1) = cmplx(0d0,0d0,kind=8)
  
      ! update B
      B = matmul(conjg(transpose(temp)),B)
      B(2,1) = cmplx(0d0,0d0,kind=8)
  
      ! update Q 
      Q = matmul(Q,temp)

    end if

  ! if B = I
  else
 
    ! compute intermediate values
    trace = A(1,1) + A(2,2)
    detm = A(1,1)*A(2,2) - A(2,1)*A(1,2)
    disc = sqrt((A(1,1)-A(2,2))**2 + 4d0*A(1,2)*A(2,1))
 
    ! compute most accurate eigenvalue
    if(abs(trace+disc) > abs(trace-disc))then
      lambda = (trace+disc)/2d0
    else
      lambda = (trace-disc)/2d0
    end if

    ! compute Q
    lambda = lambda-A(1,1)
    call z_rot3_vec4gen(dble(A(1,2)),aimag(A(1,2)),dble(lambda) &
    ,aimag(lambda),cr,ci,s,nrm)
    Q(1,1) = cmplx(cr,ci,kind=8)
    Q(2,1) = cmplx(s,0d0,kind=8)
    Q(2,2) = conjg(Q(1,1))
    Q(1,2) = -Q(2,1)

    ! update A
    A = matmul(transpose(conjg(Q)),matmul(A,Q))
    A(2,1) = cmplx(0d0,0d0,kind=8)

  end if

end subroutine z_2x2array_eig
