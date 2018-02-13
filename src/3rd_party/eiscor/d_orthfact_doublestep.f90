#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_doublestep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' doubleshift 
! algorithm on an orthogonal upper hessenberg matrix that is stored as a 
! product of givens rotations and a diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: update eigenvectors
!                    .FALSE.: no eigenvectors
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (N)
!                    array of generators for complex diagonal matrix
!
!  M               INTEGER
!                    leading dimension of Z
!
!  Z               REAL(8) array of dimension (M,N)
!                    if VEC = .TRUE. updated
!                    if VEC = .FALSE. unused
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_doublestep(VEC,N,Q,D,M,Z,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, M
  integer, intent(inout) :: ITCNT
  real(8), intent(inout) :: Q(2*(N-1)), D(N), Z(M,N)
  
  ! compute variables
  integer :: ii, ind1, ind2
  real(8) :: s1, s2
  real(8) :: b1(2), b2(2), b3(2), temp(2), nrm
  real(8) :: block(2,2), w(2,2) 

  ! compute a nonzero shift
  ! shifts +1
  if((mod(ITCNT+1,11) == 0).AND.(ITCNT.NE.0))then
    block = 0d0
    block(1,1) = 1d0 
    block(2,2) = 1d0 
    
  ! shifts -1
  else if((mod(ITCNT+1,16) == 0).AND.(ITCNT.NE.0))then
    block = 0d0
    block(1,1) = -1d0 
    block(2,2) = -1d0 
    
  ! random shift
  else if((mod(ITCNT+1,21) == 0).AND.(ITCNT.NE.0))then
    call random_number(s1)
    call random_number(s2)
    block = 0d0
    block(1,1) = s1
    block(2,2) = s2

  ! wilkinson shifts
  else

    ! get 2x2 block
    call d_orthfact_2x2diagblock(.FALSE.,Q((2*N-5):(2*N-2)),D((N-1):N),block)
      
    ! compute eigenvalues and eigenvectors
    call d_2x2array_eig(.FALSE.,block,block,w,w)
      
  end if

  ! two real shifts
  if (block(2,1).EQ.0d0) then
 
    ! build first bulge
    temp(1) = sign(1d0,block(1,1))
    temp(2) = 0d0
    call d_orthfact_buildbulge(.FALSE.,Q(1:4),D(1:2),temp,b1,b2)

    ! fusion to initialize first bulge
    b3(1) = b1(1)
    b3(2) = -b1(2)
    call d_orthfact_mergebulge(.TRUE.,Q(1:2),b3)
     
    ! first bulge update eigenvectors
    if (VEC)then
      w(1,1) = b1(1)
      w(2,1) = b1(2)
      w(1,2) = -w(2,1)
      w(2,2) = w(1,1)
      Z(:,1:2) = matmul(Z(:,1:2),w)
    end if 
    
    ! first bulge through D
    call d_rot2_swapdiag(D(1:2),b1)
      
    ! first bulge through Q
    call d_rot2_turnover(Q(1:2),Q(3:4),b1)
      
    ! second shift
    temp(1) = sign(1d0,block(2,2))
    temp(2) = 0d0
    call d_orthfact_buildbulge(.FALSE.,Q(1:4),D(1:2),temp,b2,b3)

    ! fusion to initialize second bulge
    b3(1) = b2(1)
    b3(2) = -b2(2)
    call d_orthfact_mergebulge(.TRUE.,Q(1:2),b3)
     
    ! main chasing loop
    do ii=1,(N-3)
       
      ! set ind2
      ind1 = (ii-1)
      ind2 = 2*(ii-1)
       
      ! first bulge update eigenvectors
      if (VEC)then
        w(1,1) = b1(1)
        w(2,1) = b1(2)
        w(1,2) = -w(2,1)
        w(2,2) = w(1,1)
        Z(:,(ii+1):(ii+2)) = matmul(Z(:,(ii+1):(ii+2)),w)
      end if 
        
      ! first bulge through D
      call d_rot2_swapdiag(D((ind1+2):(ind1+3)),b1)
      
      ! first bulge through Q
      call d_rot2_turnover(Q((ind2+3):(ind2+4)),Q((ind2+5):(ind2+6)),b1)
      
      ! second bulge update eigenvectors
      if (VEC)then
        w(1,1) = b2(1)
        w(2,1) = b2(2)
        w(1,2) = -w(2,1)
        w(2,2) = w(1,1)
        Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),w)
      end if 
        
      ! second bulge through D
      call d_rot2_swapdiag(D((ind1+1):(ind1+2)),b2)
      
      ! second bulge through Q
      call d_rot2_turnover(Q((ind2+1):(ind2+2)),Q((ind2+3):(ind2+4)),b2)
      
    end do
    
    ! set ind2
    ind1 = (N-3)
    ind2 = 2*(N-3)
    
    ! first bulge update eigenvectors
    if (VEC)then
      w(1,1) = b1(1)
      w(2,1) = b1(2)
      w(1,2) = -w(2,1)
      w(2,2) = w(1,1)
      Z(:,(N-1):(N)) = matmul(Z(:,(N-1):(N)),w)
    end if 
       
    ! first bulge through D
    call d_rot2_swapdiag(D((ind1+2):(ind1+3)),b1)
    
    ! first bulge fuse with Q
    call d_orthfact_mergebulge(.FALSE.,Q((ind2+3):(ind2+4)),b1)
    
    ! second bulge update eigenvectors
    if (VEC)then
      w(1,1) = b2(1)
      w(2,1) = b2(2)
      w(1,2) = -w(2,1)
      w(2,2) = w(1,1)
      Z(:,(N-2):(N-1)) = matmul(Z(:,(N-2):(N-1)),w)
    end if
        
    ! second bulge through D
    call d_rot2_swapdiag(D((ind1+1):(ind1+2)),b2)
    
    ! second bulge through Q  
    call d_rot2_turnover(Q((ind2+1):(ind2+2)),Q((ind2+3):(ind2+4)),b2)
    
    ! set ind2
    ind1 = (N-2)
    ind2 = 2*(N-2)
    
    ! last bulge update eigenvectors
    if (VEC)then
      w(1,1) = b2(1)
      w(2,1) = b2(2)
      w(1,2) = -w(2,1)
      w(2,2) = w(1,1)
      Z(:,(N-1):(N)) = matmul(Z(:,(N-1):(N)),w)
    end if 
  
    ! last bulge through D
    call d_rot2_swapdiag(D((ind1+1):(ind1+2)),b2)
    
    ! last bulge fuse with Q
    call d_orthfact_mergebulge(.FALSE.,Q((ind2+1):(ind2+2)),b2)
  
  ! complex conjugate pair
  else
  
    ! normalize shifts
    call d_rot2_vec2gen(block(1,1),sqrt(-block(2,1)*block(1,2)),temp(1),temp(2),nrm) 

    ! build bulge
    call d_orthfact_buildbulge(.TRUE.,Q(1:4),D(1:2),temp,b1,b2)

    ! turnover to initialize bulge
    temp(1) = b2(1)
    temp(2) = -b2(2)
    b3(1) = b1(1)
    b3(2) = -b1(2)
    call d_rot2_turnover(temp,b3,Q(1:2))
      
    ! fusion to finish initialization
    call d_orthfact_mergebulge(.TRUE.,Q(3:4),b3)
     
    ! update b3
    b3 = Q(1:2)
    
    ! update Q
    Q(1:2) = temp
    
    ! main chasing loop
    do ii=1,(N-3)
       
      ! set ind2
      ind1 = (ii-1)
      ind2 = 2*(ii-1)
       
      ! first bulge update eigenvectors
      if (VEC)then
        w(1,1) = b1(1)
        w(2,1) = b1(2)
        w(1,2) = -w(2,1)
        w(2,2) = w(1,1)
        Z(:,(ii+1):(ii+2)) = matmul(Z(:,(ii+1):(ii+2)),w)
      end if 
        
      ! first bulge through D
      call d_rot2_swapdiag(D((ind1+2):(ind1+3)),b1)
      
      ! first bulge through Q
      call d_rot2_turnover(Q((ind2+3):(ind2+4)),Q((ind2+5):(ind2+6)),b1)
      
      ! second bulge update eigenvectors
      if (VEC)then
        w(1,1) = b2(1)
        w(2,1) = b2(2)
        w(1,2) = -w(2,1)
        w(2,2) = w(1,1)
        Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),w)
      end if 
        
      ! second bulge through D
      call d_rot2_swapdiag(D((ind1+1):(ind1+2)),b2)
      
      ! second bulge through Q
      call d_rot2_turnover(Q((ind2+1):(ind2+2)),Q((ind2+3):(ind2+4)),b2)
      
      ! push b3 down
      call d_rot2_turnover(b3,b1,b2)
      
    ! update bulges
      temp = b2
      b2 = b3
      b3 = b1
      b1 = temp
       
    end do
    
    ! set ind2
    ind1 = (N-3)
    ind2 = 2*(N-3)
    
    ! first bulge update eigenvectors
    if (VEC)then
      w(1,1) = b1(1)
      w(2,1) = b1(2)
      w(1,2) = -w(2,1)
      w(2,2) = w(1,1)
      Z(:,(N-1):(N)) = matmul(Z(:,(N-1):(N)),w)
    end if 
       
    ! first bulge through D
    call d_rot2_swapdiag(D((ind1+2):(ind1+3)),b1)
    
    ! first bulge fuse with Q
    call d_orthfact_mergebulge(.FALSE.,Q((ind2+3):(ind2+4)),b1)
    
    ! second bulge update eigenvectors
    if (VEC)then
      w(1,1) = b2(1)
      w(2,1) = b2(2)
      w(1,2) = -w(2,1)
      w(2,2) = w(1,1)
      Z(:,(N-2):(N-1)) = matmul(Z(:,(N-2):(N-1)),w)
    end if
        
    ! second bulge through D
    call d_rot2_swapdiag(D((ind1+1):(ind1+2)),b2)
    
    ! second bulge through Q  
    call d_rot2_turnover(Q((ind2+1):(ind2+2)),Q((ind2+3):(ind2+4)),b2)
    
    ! fuse b2 and b3
    call d_orthfact_mergebulge(.FALSE.,b3,b2)
    
    ! set ind
    ind1 = (N-2)
    ind2 = 2*(N-2)
    
    ! last bulge update eigenvectors
    if (VEC)then
      w(1,1) = b3(1)
      w(2,1) = b3(2)
      w(1,2) = -w(2,1)
      w(2,2) = w(1,1)
      Z(:,(N-1):(N)) = matmul(Z(:,(N-1):(N)),w)
    end if 
  
    ! last bulge through D
    call d_rot2_swapdiag(D((ind1+1):(ind1+2)),b3)
    
    ! last bulge fuse with Q
    call d_orthfact_mergebulge(.FALSE.,Q((ind2+1):(ind2+2)),b3)
  
  end if ! end of doubleshift step

end subroutine d_orthfact_doublestep
