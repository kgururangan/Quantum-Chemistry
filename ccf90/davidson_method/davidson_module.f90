module davidson_module


        use dgeev_module, only: eig
        use dgemm_module, only: gemm, gemmt, dot
        use sort_module, only: argsort
        use print_module, only: print_matrix

        contains

                subroutine gramschmidt(Q)

                        real, intent(inout) :: Q(:,:)
                        integer :: N, i, j
                        real, allocatable :: v(:), R(:,:)

                        N = size(Q,2)
                        
                        allocate(v(N),R(N,N))

                        do j = 1,N
                           v = Q(:,j)
                           do i = 1,j-1
                              R(i,j) = dot(Q(:,i),v)
                              v = v - R(i,j)*Q(:,i)
                           end do
                           R(j,j) = norm2(v)
                           Q(:,j) = v/R(j,j)
                        end do

                        deallocate(v,R)

                end subroutine gramschmidt

                subroutine davidson(A,nroot,VR,wR,B0,tol,maxit)

                        real, intent(inout) :: A(:,:), B0(:,:)
                        real, intent(in) :: tol
                        integer, intent(in) :: nroot, maxit
                        real, allocatable, intent(out) :: VR(:,:), wR(:)
                        integer :: crsz, i, j, N, num_add, maxsz, it
                        real, allocatable :: sigma(:,:), aR(:,:), evR(:), evI(:), G(:,:),&
                                             rv(:), v(:), q(:), B(:,:), addB(:,:)
                        integer, allocatable :: id(:)
                        real :: resid, s, resv(nroot)
                        integer, parameter :: nvec = 20
                        real, parameter :: threshold = 1e-05

                        ! matrix dimension
                        N = size(A,2)

                        ! current subspace size
                        crsz = size(B0,2)

                        ! maximum subspace size dimension
                        maxsz = nvec * nroot

                        ! allocate the output right eigenvectors and eigenvalues
                        allocate(VR(N,nroot),wR(nroot),rv(N),q(N),B(N,maxsz),addB(N,nroot),v(N))

                        ! populate the B matrix with the initial guess
                        do i = 1,crsz
                           B(:,i) = B0(:,i)
                        end do

                        ! set up the residuals to be large upon entering the Davidson loop,
                        ! for the sake of the first iteration
                        do i = 1,nroot
                           resv(i) = 100.0
                        end do

                        ! main davidson update loop
                        do it = 1,maxit

                                write(*,fmt=*) ''
                                write(*,fmt=*) 'Iteration - ',it,'Subspace dim = ',crsz
                                write(*,fmt=*) '------------------------------------------------------------------------'

                                ! allocate arrays
                                allocate(sigma(N,crsz),&
                                         aR(crsz,crsz),evR(crsz),evI(crsZ),G(crsz,crsz),&
                                         id(crsz))

                                ! orthogonalize guess space
                                call gramschmidt(B(:,1:crsz))

                                ! matrix-vector product
                                call gemm(A,B(:,1:crsz),sigma)

                                ! create projected subspace interaction matrix
                                !call gemm(transpose(B(:,1:crsz)),sigma,G)
                                call gemmt(B(:,1:crsz),sigma,G)

                                ! diagonalize projected eigenvalue problem
                                call eig(G,aR,evR,evI)

                                ! sort eigenvalues
                                id = argsort(evR)
                                evR = evR(id)
                                aR = aR(:,id)

                                resid = 0.0
                                num_add = 0
                                do j = 1,nroot

                                   ! if root j is already converged, don't touch it
                                   if (resv(j) < tol) then
                                      cycle
                                   else

                                      ! calculate ritz vector
                                      v = 0.0
                                      rv = 0.0
                                      do i = 1,crsz
                                         v = v + B(:,i)*aR(i,j)
                                         rv = rv + sigma(:,i)*aR(i,j)
                                      end do

                                      ! store current eigenvector/value pair
                                      VR(:,j) = v
                                      wR(j) = evR(j)

                                      ! calculate residual vector
                                      rv = rv - v*evR(j)

                                      ! increment residuum
                                      resv(j) = norm2(rv)
                                      resid = resid + resv(j)

                                      ! calculate DPR updated vector
                                      q = rv/(evR(j)-A(j,j))
                                      q = q/norm2(q)

                                      ! orthogonalize against all vectors in current search space
                                      do i = 1,crsz
                                         v = B(:,i)/norm2(B(:,i))
                                         s = dot(v,q)
                                         q = q - s*v
                                      end do
                                      ! and those that are to be added to the search space
                                      do i = 1,num_add
                                         v = addB(:,i)/norm2(addB(:,i))
                                         s = dot(v,q)
                                         q = q - s*v
                                      end do

                                      ! add q to search space if its norm is large enough
                                      if (norm2(q) > threshold) then
                                          num_add = num_add + 1
                                          addB(:,num_add) = q
                                      end if

                                   end if

                                end do

                                ! print status of each root
                                do j = 1,nroot
                                        write(*,fmt=*) 'Root - ',j,'e = ',wR(j),'|r| = ',resv(j)
                                end do

                                ! check exit condition
                                if (resid <= tol) then
                                   write(*,fmt='(/a)') 'Davidson converged!'
                                   do i = 1,nroot
                                      print*,'Root',i,'=',wR(i)
                                   end do
                                   exit

                                ! otherwise, check if number of search vectors exceeds max size
                                elseif (crsz+num_add > maxsz) then

                                   ! collapse subspace    
                                   call gemm(B(:,1:crsz),aR(1:crsz,1:nroot),addB)
                                   do i = 1,nroot
                                      B(:,i) = addB(:,i)
                                   end do

                                   ! update search subspace dimension
                                   crsz = nroot

                                ! if crsz+addB <= maxsz, expand the subspace
                                else 

                                    do i = 1,num_add
                                       B(:,crsz+i) = addB(:,i)
                                    end do

                                    ! update search subspace dimension
                                    crsz = crsz + num_add

                                end if

                                ! deallocate temporary arrays
                                deallocate(sigma,aR,evR,evI,G,id)

                        end do

                        ! deallocate rest of the arrays
                        deallocate(rv,q,B,addB,v)

                end subroutine davidson

                !subroutine update()

                !end subroutine update()


end module davidson_module
