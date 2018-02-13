#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_2x2diagblocks
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a set of two by two diagonal blocks of a 
! unitary plus rank one matrix pencil stored in factored form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  TOP             LOGICAL
!                    .TRUE.: top block is computed
!                    .FALSE.: bottom block is computed
!
!  HESS            LOGICAL
!                    .TRUE.: the extended hessenberg part is included
!                    .FALSE.: only upper triangular parts are returned
!
!  QZ              LOGICAL
!                    .TRUE.: second triangular factor is assumed nonzero
!                    .FALSE.: second triangular factor is assumed to be identity
!
!  P               LOGICAL
!                    position flag for Q
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for first sequence of rotations
!                    if HESS = .FALSE., unused
!
!  D1,D2           REAL(8) arrays of dimension (4)
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!                    if QZ = .FALSE., D2 is unused
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (6)
!                    array of generators for upper-triangular parts of the pencil
!                    if QZ = .FALSE., C2 and B2 are unused
!
! OUTPUT VARIABLES:
!
!  A,B             COMPLEX(8) array of dimension (2,2)
!                    on exit contains the desired 2x2 block from
!                    the extended hessenberg matrix 
!                    if QZ = .FALSE., B is unused
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_2x2diagblocks(TOP,HESS,QZ,P,Q,D1,C1,B1,D2,C2,B2,A,B)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: TOP, HESS, QZ, P
  real(8), intent(in) :: Q(6), D1(4), C1(6), B1(6)
  real(8), intent(in) :: D2(4), C2(6), B2(6)
  complex(8), intent(inout) :: A(2,2), B(2,2)
  
  ! compute variables
  complex(8) :: H(2,2)
  
  ! compute A
  A = cmplx(0d0,0d0,kind=8)

  ! first column of T
  A(1,1) = cmplx(-B1(3)/C1(3),0d0,kind=8)

  ! second column of T
  A(2,2) = cmplx(-B1(6)/C1(6),0d0,kind=8)
  A(1,2) = (cmplx(-B1(1),B1(2),kind=8)*cmplx(B1(4),B1(5),kind=8) &
      + A(2,2)*cmplx(C1(1),C1(2),kind=8)*cmplx(C1(4),-C1(5),kind=8))/cmplx(C1(3),0d0,kind=8)
  
  ! apply diagonal
  A(1,:) = cmplx(D1(1),D1(2),kind=8)*A(1,:)
  A(2,:) = cmplx(D1(3),D1(4),kind=8)*A(2,:)
  
  ! extended hessenberg part
  if (HESS.AND.TOP) then

    ! build local Q
    H(1,1) = cmplx(Q(1),Q(2),kind=8)
    H(2,1) = cmplx(Q(3),0d0,kind=8)
    H(1,2) = cmplx(-Q(3),0d0,kind=8)
    H(2,2) = cmplx(Q(1),-Q(2),kind=8)
    
    ! include adjacent rotation
    if (P) then  
      H(2,:) = cmplx(Q(4),Q(5),kind=8)*H(2,:)
    else
      H(:,2) = cmplx(Q(4),Q(5),kind=8)*H(:,2)
    end if
    
    ! set output
    A = matmul(H,A)
    
  else if (HESS) then

    ! build local Q
    H(1,1) = cmplx(Q(4),Q(5),kind=8)
    H(2,1) = cmplx(Q(6),0d0,kind=8)
    H(1,2) = cmplx(-Q(6),0d0,kind=8)
    H(2,2) = cmplx(Q(4),-Q(5),kind=8)
    
    ! include adjacent rotation
    if (P) then  
      H(:,1) = cmplx(Q(1),Q(2),kind=8)*H(:,1)
    else
      H(1,:) = cmplx(Q(1),Q(2),kind=8)*H(1,:)
    end if
    
    ! set output
    A = matmul(H,A)
    
  end if

  ! compute B
  if (QZ) then
  
    ! initialize 
    B = cmplx(0d0,0d0,kind=8)

    ! first column of T
    B(1,1) = cmplx(-B2(3)/C2(3),0d0,kind=8)
        
    ! second column of T
    B(2,2) = cmplx(-B2(6)/C2(6),0d0,kind=8)
    B(1,2) = (cmplx(-B2(1),B2(2),kind=8)*cmplx(B2(4),B2(5),kind=8) &
        + B(2,2)*cmplx(C2(1),C2(2),kind=8)*cmplx(C2(4),-C2(5),kind=8))/cmplx(C2(3),0d0,kind=8)
    
    ! apply diagonal
    B(1,:) = cmplx(D2(1),D2(2),kind=8)*B(1,:)
    B(2,:) = cmplx(D2(3),D2(4),kind=8)*B(2,:)  
    
  end if

end subroutine z_upr1fact_2x2diagblocks
