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

                subroutine gemmt(A,B,C)

                       real, intent(in) :: A(:,:), B(:,:)
                       real, intent(out) :: C(:,:)
                       integer :: m, n, k
                       real :: alpha = 1.0, beta = 0.0
                       m = ubound(A,2)
                       n = ubound(B,2)
                       k = ubound(A,1)
                       call dgemm('t','n',m,n,k,alpha,A,k,B,k,beta,C,m)

                end subroutine gemmt

                subroutine gemv(A,b,c)

                        real, intent(in) :: A(:,:), b(:)
                        real, intent(out) :: c(:)
                        real :: alpha = 1.0, beta = 0.0
                        integer :: m, n

                        m = size(A,1)
                        n = size(A,2)

                        call dgemv('n',m,n,alpha,A,m,b,1,beta,c,1)

                end subroutine gemv

                function dot(a,b) result(x)

                        real, intent(in) :: a(:), b(:)
                        real :: x
                        integer :: N, i

                        N = size(a)

                        x = 0.0
                        do i = 1,N
                           x = x + a(i) * b(i)
                        end do

                end function dot


end module dgemm_module
