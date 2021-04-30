module diis

      use dgemm_module, only: dot

      implicit none

      contains

                subroutine diis_xtrap(X_xtrap, X_list, diis_resid_list)

                        real, intent(in) :: X_list(:,:), diis_resid_list(:,:)
                        real, intent(inout) :: X_xtrap(size(X_list,1))
                        real :: B(size(X_list,2)+1,size(X_list,2)+1), rhs(size(X_list,2)+1), coeff(size(X_list,2)+1)
                        integer :: i, j, k, diis_size, vec_dim, B_dim

                        vec_dim = size(X_list,1)
                        diis_size = size(X_list,2)
                        B_dim = diis_size + 1

                        B = 0.0
                        B(B_dim,:) = -1.0
                        B(:,B_dim) = -1.0
                        B(B_dim,B_dim) = 0.0

                        do i = 1,diis_size
                           do j = i,diis_size
                              B(i,j) = dot(diis_resid_list(:,i),diis_resid_list(:,j))
                              B(j,i) = B(i,j)
                           end do
                        end do

                        rhs = 0.0
                        rhs(B_dim) = -1.0
                        call solvegauss(B,rhs,B_dim,coeff)
                        X_xtrap = 0.0
                        do k = 1,diis_size
                           X_xtrap = X_xtrap + coeff(k)*X_list(:,k)
                        end do
                end subroutine diis_xtrap                                        


                subroutine solvegauss(A,b,N,x)

                        integer, intent(in) :: N
                        real, intent(inout) :: A(N,N), b(N)
                        real, intent(out) :: x(N)
                        real :: m, xsum
                        integer :: i, j

                        do j = 1,N-1
                           do i = N,j+1,-1
                              m = A(i,j)/A(j,j)
                              A(i,:) = A(i,:) - m*A(j,:)
                              b(i) = b(i) - m*b(j)
                           end do
                        end do

                        x = 0.0
                        x(N) = b(N)/A(N,N)
                        do i = N-1,1,-1
                           xsum = 0.0
                           do j = N,i+1,-1
                              xsum = xsum + A(i,j)*x(j)
                           end do
                           x(i) = (b(i)-xsum)/(A(i,i))
                        end do

                end subroutine solvegauss



end module diis
