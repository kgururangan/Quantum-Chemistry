module davidson_module


        use dgeev_module, only: eig
        use dgemm_module, only: gemm, gemmt, gemv, dot
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
                        integer :: crsz, i, j, k, N, num_add, maxsz, it
                        real, allocatable :: sigma(:,:), aR(:,:), evR(:), evI(:), G(:,:),&
                                             rv(:), v(:), q(:), B(:,:), addB(:,:)
                        integer, allocatable :: id(:)
                        real :: resid, s, resv(nroot)
                        integer, parameter :: nvec = 20
                        real, parameter :: threshold = 1.0e-05

                        ! matrix dimension
                        N = size(A,2)

                        ! current subspace size
                        crsz = size(B0,2)

                        ! maximum subspace size dimension
                        maxsz = nvec * nroot

                        ! allocate the output right eigenvectors and eigenvalues
                        allocate(VR(N,nroot),wR(nroot),rv(N),q(N),B(N,maxsz),addB(N,nroot),sigma(N,maxsz),v(N))

                        ! populate the B matrix with the initial guess
                        B(:,1:crsz) = B0

                        ! set up the residuals to be large upon entering the Davidson loop,
                        ! for the sake of the first iteration
                        resv(1:nroot) = 100.0
                        num_add = nroot

                        ! main davidson update loop
                        do it = 1,maxit

                                write(*,fmt=*) ''
                                write(*,'(1x,a13,1x,i0,10x,a15,1x,i0)') 'Iteration - ',it,'Subspace dim - ',crsz
                                write(*,fmt=*) '-----------------------------------------------------------'
                                
                                ! allocate arrays
                                allocate(aR(crsz,crsz),evR(crsz),evI(crsZ),G(crsz,crsz),id(crsz))

                                ! orthogonalize guess space
                                call gramschmidt(B(:,1:crsz))

                                ! matrix-vector product
                                k = max(crsz-num_add+1,1)
                                call gemm(A,B(:,k:crsz),sigma(:,k:crsz))

                                ! create projected subspace interaction matrix
                                call gemmt(B(:,1:crsz),sigma(:,1:crsz),G)

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
                                      call gemv(B(:,1:crsz),aR(:,j),v)
                                      call gemv(sigma(:,1:crsz),aR(:,j),rv)

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
                                   write(*,'(1x,a8,i0,8x,a4,f12.8,8x,a6,f12.8)') 'Root - ',j,'e = ',wR(j),'|r| = ',resv(j)
                                end do

                                ! check exit condition
                                if (resid <= tol) then
                                   write(*,fmt='(/a)') 'Davidson converged!'
                                   do i = 1,nroot
                                      write(*,'(1x,a4,2x,i0,a4,f12.8)') 'Root',i,' = ',wR(i)
                                   end do
                                   exit

                                ! otherwise, check if number of search vectors exceeds max size
                                elseif (crsz+num_add > maxsz) then

                                   ! collapse subspace    
                                   call gemm(B(:,1:crsz),aR(1:crsz,1:nroot),addB)
                                   B(:,1:nroot) = addB
                                   ! update search subspace dimension
                                   crsz = nroot

                                ! if crsz+addB <= maxsz, expand the subspace
                                else 

                                    B(:,crsz+1:crsz+num_add) = addB(:,1:num_add)
                                    ! update search subspace dimension
                                    crsz = crsz + num_add

                                end if

                                ! deallocate temporary arrays
                                deallocate(aR,evR,evI,G,id)

                        end do

                        ! deallocate rest of the arrays
                        deallocate(rv,q,B,sigma,addB,v)

                end subroutine davidson


                subroutine calc_left(A,VR,wR,tol,maxit,shift,VL)

                        use diis, only: diis_xtrap

                        real, intent(in) :: A(:,:), VR(:,:), wR(:), tol, shift
                        integer, intent(in) :: maxit
                        real, allocatable, intent(out) :: VL(:,:)
                        integer, parameter :: ndiis = 6
                        integer :: N, nroot, it, it_macro, i, j
                        real, allocatable :: L(:), L_list(:,:), L_resid_list(:,:), rv(:), LA(:), wL(:)
                        real :: denom, res, LR

                        N = size(A,1)
                        nroot = size(VR,2)

                        allocate(VL(N,nroot),L(N),LA(N),L_list(N,ndiis),L_resid_list(N,ndiis),rv(N),wL(nroot))
                        
                        do i = 1,nroot
                                
                                write(*,fmt=*) ''
                                write(*,fmt=*) 'Jacobi solve for root - ',i
                        write(*,fmt=*) '               Iteration      E(left)                       Residuum'         
                        write(*,fmt=*) '------------------------------------------------------------------------'
                                L = vR(:,i)
                                L_list = 0.0
                                L_resid_list = 0.0
                                it_macro = 0

                                do it = 0,maxit
                                   
                                   call gemv(transpose(A),-1.0*L,rv)
                                   do j = 1,N
                                      denom = A(j,j)-wR(i)+shift
                                      L(j) = (A(j,j)+shift)*L(j)/denom + rv(j)/denom
                                   end do
                                   call gemv(transpose(A),L,LA)

                                   rv = LA - wR(i)*L
                                   res = norm2(rV)

                                   wL(i) = norm2(LA)/norm2(L)

                                   L_list(:,mod(it,ndiis)+1) = L
                                   L_resid_list(:,mod(it,ndiis)+1) = rv

                                   if (mod(it+1,ndiis) == 0) then
                                      it_macro = it_macro + 1
                                      write(*,fmt='(1x,a15,1x,i0)') 'DIIS cycle - ',it_macro
                                      call diis_xtrap(L,L_list,L_resid_list)
                                   end if
                                    
                                   LR = dot(L,VR(:,i))
                                   L = L/LR

                                   write(*,fmt=*) it,' ',wL(i),' ',res

                                   if (res < tol) then
                                      write(*,fmt=*) 'Root',i,'converged!'
                                      exit
                                   end if

                                end do

                           end do

                   end subroutine calc_left



end module davidson_module
