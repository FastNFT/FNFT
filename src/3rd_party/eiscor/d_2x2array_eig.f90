#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_2x2array_eig 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the real Schur decomposition of a general
! 2x2 double matrix pencil.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies B is not the identity
!                    .FALSE. implies B is the identity matrix
!
!  A,B             REAL(8) array of dimension (2,2)
!                    Contains the 2x2 matrix. On exit (A,B) are real and
!                    quasi-uppertriangular. If FLAG=.FALSE. B is
!                    unused.
!
! OUTPUT VARIABLES:
!
!  Q,Z             COMPLEX(8) array of dimension (2,2)
!                    On exit the columns of Q and Z contain the left and right
!                    Schur vectors of the pencil (A,B). If FLAG=.FALSE. Z is
!                    unused.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_2x2array_eig(FLAG,A,B,Q,Z)
  
  implicit none
  
  ! input variables
  logical, intent(inout) :: FLAG
  real(8), intent(inout) :: A(2,2), B(2,2), Q(2,2), Z(2,2)
  
  ! compute variables
  real(8) :: trace, disc, detm, temp
  real(8) :: c, s, nrm, T(2,2)
 
  ! B not identity
  if (FLAG) then
  
    ! make B upper-triangular
    ! create eliminator
    call d_rot2_vec2gen(B(1,1),B(2,1),c,s,nrm)
      
    ! store in Q
    Q(1,1) = c
    Q(2,1) = s
    Q(2,2) = Q(1,1)
    Q(1,2) = -Q(2,1)
      
    ! update B
    B = matmul(transpose(Q),B)
    B(2,1) = 0d0
      
    ! update A
    A = matmul(transpose(Q),A)
      
    ! set Z
    Z(1,1) = 1d0
    Z(2,1) = 0d0
    Z(2,2) = 1d0
    Z(1,2) = 0d0
    
    ! check for infinite eigenvalues
    if (abs(B(1,1)).EQ.0d0) then
  
      ! create eliminator
      call d_rot2_vec2gen(A(1,1),A(2,1),c,s,nrm)
      
      ! store in T
      T(1,1) = c
      T(2,1) = s
      T(2,2) = T(1,1)
      T(1,2) = -T(2,1)
      
      ! update A
      A = matmul(transpose(T),A)
      A(2,1) = 0d0
      
      ! update B
      B = matmul(transpose(Q),B)
      
      ! update Q
      Q = matmul(Q,T)
      
      ! return
      return
  
    ! check for infinite eigenvalues
    else if (abs(B(2,2)).EQ.0d0) then
  
      ! create eliminator
      call d_rot2_vec2gen(A(2,2),A(2,1),c,s,nrm)
      
      ! store in T
      T(1,1) = c
      T(2,1) = s
      T(2,2) = T(1,1)
      T(1,2) = -T(2,1)
      
      ! update A
      A = matmul(A,T)
      A(2,1) = 0d0
      
      ! update B
      B = matmul(B,T)
      
      ! update Z
      Z = matmul(Z,T)
      
      ! return
      return

    ! no infinite eigenvalues
    else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this section needs finishing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end if 
  ! B identity
  else  

    ! compute intermediate values
    trace = A(1,1) + A(2,2)
    detm = A(1,1)*A(2,2) - A(2,1)*A(1,2)
    temp = (A(1,1)-A(2,2))**2 + 4d0*A(1,2)*A(2,1)
    
    ! imaginary roots
    if (temp.LT.0d0) then
      
      ! move A to standard form (A(1,1) = A(2,2))
   
      ! compute Q
      temp = (A(2,2)-A(1,1))
      if (temp.NE.0d0) then
        
        temp = (A(1,2)+A(2,1))/temp

        if (abs(temp).GE.1d0) then
          temp = temp*(1d0+sqrt(1d0+1d0/temp/temp))
          call d_rot2_vec2gen(temp,1d0,Q(1,1),Q(2,1),temp)
        else 
          temp = temp+sign(1d0,temp)*sqrt(1d0+temp*temp)
          call d_rot2_vec2gen(temp,1d0,Q(1,1),Q(2,1),temp)
        end if

      else

        Q(1,1) = 1d0
        Q(2,1) = 0d0

      end if

      Q(2,2) = Q(1,1)
      Q(1,2) = -Q(2,1)

      ! update A
      A = matmul(A,Q)
      A = matmul(transpose(Q),A)

    ! real roots
    else
     
      ! sqrt of discriminant 
      disc = sqrt(temp)
      
      ! compute most accurate eigenvalue
      if(abs(trace+disc) > abs(trace-disc))then
        temp = (trace+disc)/2d0
      else
        temp = (trace-disc)/2d0
      end if

      ! compute Schur vectors
      call d_rot2_vec2gen(A(1,2),temp-A(1,1),Q(1,1),Q(2,1),trace)
      Q(2,2) = Q(1,1)
      Q(1,2) = -Q(2,1)
      
      ! update A
      A = matmul(transpose(Q),matmul(A,Q))
      A(2,1) = 0d0
      
    end if
  
  end if
  
end subroutine d_2x2array_eig
