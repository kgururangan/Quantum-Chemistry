module dgemm_module


        implicit none

        contains

                subroutine gemm(A,B,C)

                       real, intent(in) :: A(:,:), B(:,:)
                       real, intent(out) :: C(:,:)
                       integer :: m, n, k
                       real :: alpha = 1.0, beta = 0.0
                       m = ubound(A,1)
                       n = ubound(B,2)
                       k = ubound(B,1)
                       call dgemm('n','n',m,n,k,alpha,A,m,B,k,beta,C,m)

                end subroutine gemm

end module dgemm_module
