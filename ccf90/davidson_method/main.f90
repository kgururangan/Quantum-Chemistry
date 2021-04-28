program main

        use davidson_module, only: davidson
        use print_module, only: print_matrix
        use dgemm_module, only: gemm, gemmt

        implicit none

        integer, parameter :: N = 1000, nroot = 5, maxit=3000
        real, parameter :: sparsity = 1e-01, herm = 1, tol=1e-06
        real :: A(N,N)
        real, allocatable :: B0(:,:), VR(:,:), wR(:)

        call get_test_matrix(N,sparsity,herm,A)

        call get_diagonal_guess(A,nroot,B0)

        call davidson(A,nroot,VR,wR,B0,tol,maxit)

        ! Do this last since A changes upon output!
        write(*,'(/a)') 'Exact answer : '
        call dense_diag(A,nroot)


        contains

                subroutine get_test_matrix(N,sparsity,herm,A)

                        integer, intent(in) :: N
                        real, intent(in) :: sparsity, herm
                        real, intent(out) :: A(N,N)
                        real, parameter :: lb = -2.0, ub = 2.0
                        integer :: i, j
                        real :: r

                        do i = 1,N
                           do j = 1,N
                              call random_number(r)
                              A(i,j) = sparsity*((ub-lb)*r + lb)
                           end do
                        end do

                        do i = 1,N
                           A(i,i) = A(i,i) + i
                        end do

                        A = herm*A + (1.0-herm)*transpose(A)

                 end subroutine get_test_matrix

                 subroutine dense_diag(A,nroot)

                        use dgeev_module, only: eig
                        use sort_module, only: argsort

                        real, intent(in) :: A(:,:)
                        integer, intent(in) :: nroot
                        integer :: N, i
                        integer, allocatable :: idx(:)
                        real, allocatable :: VR(:,:), VL(:,:), wR(:), wI(:)

                        N = size(A,1)
                        allocate(VR(N,N), VL(N,N), wR(N), wI(N), idx(N))

                        call eig(A,VR,wR,wI)
                        idx = argsort(wR)
                        wR = wR(idx)
                        do i = 1,nroot
                           print*,'Eigenvalue',i,'=',wR(i)
                        end do
                        deallocate(VR,VL,wR,wI,idx)

                 end subroutine dense_diag
        
                         
                 subroutine get_diagonal_guess(A,nroot,B0)

                         use sort_module, only: argsort
        
                         real, intent(in) :: A(:,:)
                         integer, intent(in) :: nroot
                         real, allocatable :: B0(:,:), D(:)
                         integer :: N, i
                         integer, allocatable :: idx(:)

                         N = size(A,2)

                         allocate(B0(N,nroot),D(N),idx(N))

                         do i = 1,N
                            D(i) = A(i,i)
                         end do

                         idx = argsort(D)

                         B0 = 0.0 * B0
                         do i = 1,nroot
                            B0(idx(i),i) = 1.0
                         end do

                         deallocate(D,idx)

                 end subroutine get_diagonal_guess


end program main
                        
                         


