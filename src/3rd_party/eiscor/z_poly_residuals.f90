!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_poly_residuals 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the roots of a polynomial expressed in the 
! monomial basis using the fast algorithm described in:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!  ROOTS           COMPLEX(8) array of dimension (N)
!                    computed roots
!
! OUTPUT VARIABLES:
!
!  RESIDUALS       REAL(8) array of dimension (N)
!                    residuals
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_poly_residuals(N,COEFFS,ROOTS,K,RESIDUALS)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: K
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)
  real(8), intent(inout) :: RESIDUALS(N)
  
  ! compute variables
  integer :: ii, jj, kk
  complex(8) :: p, dp 

  ! check K
  if (K<0) then
    K = 0
  end if

  ! newton corrections
  do kk = 1,K

    ! loop through roots
    do ii = 1,N

      ! inside unit circle
      if (abs(ROOTS(ii)) <= 1d0) then

        ! horners rule
        p = COEFFS(1)
        dp = dble(N)*p      
        do jj = 1,(N-1)
          p = ROOTS(ii)*p + COEFFS(jj+1)
          dp = ROOTS(ii)*dp + dble(N-jj)*COEFFS(jj+1)
        end do
        p = ROOTS(ii)*p + COEFFS(N+1)

        ! store residual	
        ROOTS(ii) = ROOTS(ii) - p/dp

      ! outside the unit circle 
      else

        ! horners rule
        p = COEFFS(N+1)
        dp = COEFFS(N)
        do jj = 1,(N-1)
          p = p/ROOTS(ii) + COEFFS(N+1-jj)
          dp = dp/ROOTS(ii) + dble(jj+1)*COEFFS(N-jj)
        end do
        p = p/ROOTS(ii) + COEFFS(1)

        ! store residual	
        ROOTS(ii) = ROOTS(ii) - ROOTS(ii)*p/dp

      end if

    end do 

  end do

  ! loop through roots
  do ii = 1,N

    ! inside unit circle
    if (abs(ROOTS(ii)) <= 1d0) then

      ! horners rule
      p = COEFFS(1)
      dp = dble(N)*p      
      do jj = 1,(N-1)
        p = ROOTS(ii)*p + COEFFS(jj+1)
        dp = ROOTS(ii)*dp + dble(N-jj)*COEFFS(jj+1)
      end do
      p = ROOTS(ii)*p + COEFFS(N+1)

      ! store residual	
      RESIDUALS(ii) = abs(p/dp)

    ! outside the unit circle 
    else

      ! horners rule
      p = COEFFS(N+1)
      dp = COEFFS(N)
      do jj = 1,(N-1)
        p = p/ROOTS(ii) + COEFFS(N+1-jj)
        dp = dp/ROOTS(ii) + dble(jj+1)*COEFFS(N-jj)
      end do
      p = p/ROOTS(ii) + COEFFS(1)

      ! store residual	
      RESIDUALS(ii) = abs(ROOTS(ii)*p/dp)

    end if

  end do 

end subroutine z_poly_residuals
