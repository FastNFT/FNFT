!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_decompress
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine decompresses a factored unitary plus rank one 
! (upr1fact) matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for the unitary rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  C               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of rotations
!
!  B               REAL(8) array of dimension (3*N)
!                    array of generators for second sequence of rotations
!
! OUTPUT VARIABLES:
!
!  H               COMPLEX(8) array of dimension (N,N)
!                    extended hessenberg matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_decompress(N,P,Q,D,C,B,H)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  real(8), intent(in) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  complex(8), intent(inout) :: H(N,N)
  
  ! compute variables
  integer :: ii, jj, ind, rowind
  complex(8) :: temp(2,2)
  
  ! initialize H
  H = cmplx(0d0,0d0,kind=8)

  ! decompress triangular part
  call z_upr1utri_decompress(.FALSE.,N,D,C,B,H)

  ! loop through position flags   
  ind = 1
  do ii = 1,N-2

    ! if P(ii) == .TRUE. then we work 
    ! backwards to multiply in the corresponding Qs
    if ( P(ii) ) then

      ! loop backwards through Q factors
      do jj = 1,ind

        ! set rowind
        rowind = ii + 1 - jj
        
        ! multiply in Q factor
        temp(1,1) = cmplx(Q(3*rowind-2),Q(3*rowind-1),kind=8)
        temp(2,1) = cmplx(Q(3*rowind),0d0,kind=8)
        temp(1,2) = -temp(2,1)
        temp(2,2) = conjg(temp(1,1))
        H(rowind:rowind+1,:) = matmul(temp,H(rowind:rowind+1,:))

        ! reset ind
        ind = 1

      end do

    ! if P(ii) == .FALSE. then we increment ind
    else
   
      ind = ind + 1

    end if 

  end do

  ! final do loop
  do jj = 1,ind

    ! set rowind
    rowind = N - jj
    
    ! multiply in Q factor
    temp(1,1) = cmplx(Q(3*rowind-2),Q(3*rowind-1),kind=8)
    temp(2,1) = cmplx(Q(3*rowind),0d0,kind=8)
    temp(1,2) = -temp(2,1)
    temp(2,2) = conjg(temp(1,1))
    H(rowind:rowind+1,:) = matmul(temp,H(rowind:rowind+1,:))

  end do

end subroutine z_upr1fact_decompress
