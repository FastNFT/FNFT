subroutine my_matmul(M, N, A, B)
      implicit none
      integer, intent(in) :: M, N
      complex(8), intent(in) :: A(M,M)
      complex(8), intent(inout) :: B(M,N)
      B = matmul(A, B)
end subroutine my_matmul
