module davidson_module


        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use dgeev_module, only: eig
        use dgemm_module, only: gemm, gemmt, dot
        use sort_module, only: argsort

        contains

                subroutine gramschmidt(Q)
                ! Description: orthogonalizes the columns of input 2D array Q

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

                subroutine davidson_cis(sys,fA,fB,vA,vB,vC,N,nroot,VR,wR,tol,maxit)

                        use cis_updates, only: HC_matmat

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        integer, intent(in) :: N
                        real, intent(in) :: tol
                        integer, intent(in) :: nroot, maxit
                        real, intent(out) :: VR(N,nroot), wR(nroot)
                        integer :: crsz, i, j, a, k, num_add, maxsz, it, cnt, idx(N), crsz0
                        real, allocatable :: sigma(:,:), aR(:,:), evR(:), evI(:), G(:,:),&
                                             rv(:), v(:), q(:), B(:,:), addB(:,:)
                        integer, allocatable :: id(:)
                        real :: resid, s, resv(nroot), D1(N), B0(N,nroot)
                        integer, parameter :: nvec = 10
                        real, parameter :: threshold = 1.0e-05, hatoev = 27.2114

                        ! print header
                        write(*,fmt=*) 'CIS ALGORITHM - BLOCK DAVIDSON SOLVER'
                        write(*,fmt=*) 'NUMBER OF ROOTS - ',nroot
                        write(*,fmt=*) 'MAX SUBSPACE SZ - ',nroot*nvec
                        write(*,fmt=*) 'CONVRG THRESH = ',tol
                        write(*,fmt=*) 'SUBSPACE EXPAND THRESH = ',threshold

                        ! get the diagonal guess
                        cnt = 1
                        do i = 1,sys%Nocc_a
                           do a = 1,sys%Nunocc_a
                              D1(cnt) = fA%uu(a,a)-fA%oo(i,i)+vA%ouuo(i,a,a,i)
                              cnt = cnt + 1
                           end do
                        end do
                        do i = 1,sys%Nocc_b
                           do a = 1,sys%Nunocc_b
                              D1(cnt) = fB%uu(a,a)-fB%oo(i,i)+vC%ouuo(i,a,a,i)
                              cnt = cnt + 1
                           end do
                        end do

                        idx = argsort(D1)
                        B0 = 0.0
                        do i = 1,nroot
                           B0(idx(i),i) = 1.0
                        end do

                        ! get current subspace size
                        crsz = size(B0,2)
                        ! initialize previous subspace size =/ crsz at first iteration
                        crsz0 = 0

                        ! maximum subspace size dimension
                        maxsz = nvec * nroot

                        ! allocate the output right eigenvectors and eigenvalues
                        allocate(rv(N),q(N),B(N,maxsz),sigma(N,maxsz),addB(N,nroot),v(N))

                        ! populate the B matrix with the initial guess
                        B(:,1:crsz) = B0

                        ! For the sake of the first iteration:
                        ! (1) initialize residuals to be large upon entering the Davidson loop
                        resv(1:nroot) = 100.0
                        ! (2) set num_add equal to nroot
                        num_add = nroot

                        ! main davidson loop
                        do it = 1,maxit

                                write(*,fmt=*) ''
                                write(*,'(1x,a13,1x,i0,10x,a15,1x,i0)') 'Iteration - ',it,'Subspace dim - ',crsz
                                write(*,fmt=*) '-----------------------------------------------------------'

                                ! allocate arrays
                                allocate(aR(crsz,crsz),evR(crsz),evI(crsZ),G(crsz,crsz),id(crsz))

                                ! check if subspace is stangated
                                if (crsz0 == crsz) then
                                   write(*,fmt=*) 'Davidson stagnated!'
                                   exit
                                else
                                   crsz0 = crsz
                                end if

                                ! orthogonalize guess space
                                call gramschmidt(B(:,1:crsz))

                                ! matrix-vector product
                                k = max(crsz-num_add+1,1)
                                call HC_matmat(sys,fA,fB,vA,vB,vC,B(:,k:crsz),N,crsz-k+1,sigma(:,k:crsz))

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
                                      call update_C(rv,sys,fA,fB,vA,vB,vC,evR(j),q)

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
                                      write(*,'(1x,a4,2x,i0,a4,f12.8,2x,a3)') 'Root',i,' = ',wR(i)*hatoev,'eV'
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
                                deallocate(aR,evR,evI,G,id)

                        end do

                        ! deallocate rest of the arrays
                        deallocate(rv,q,B,sigma,addB,v)

                end subroutine davidson_cis

                subroutine update_C(q,sys,fA,fB,vA,vB,vC,omega,qout)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: q(:), omega
                        real, intent(out) :: qout(sys%Nocc_a*sys%Nunocc_a+sys%Nocc_b*sys%Nunocc_b)
                        real :: c1a(sys%Nunocc_a,sys%Nocc_a), c1b(sys%Nunocc_b,sys%Nocc_b)
                        integer :: a, i

                        c1a = reshape(q(1:sys%Nocc_a*sys%Nunocc_a),shape(c1a))
                        c1b = reshape(q(sys%Nocc_a*sys%Nunocc_a+1:sys%Nocc_a*sys%Nunocc_a+sys%Nocc_b*sys%Nunocc_b),shape(c1b))

                        do i = 1,sys%Nocc_a
                           do a = 1,sys%Nunocc_a
                              denom = fA%uu(a,a)-fA%oo(i,i)+vA%ouuo(i,a,a,i)
                              c1a(a,i) = c1a(a,i)/(omega-denom)
                           end do
                        end do

                        do i = 1,sys%Nocc_b
                           do a = 1,sys%Nunocc_b
                              denom = fB%uu(a,a)-fB%oo(i,i)+vC%ouuo(i,a,a,i)
                              c1b(a,i) = c1b(a,i)/(omega-denom)
                           end do
                        end do

                        qout(1:sys%Nocc_a*sys%Nunocc_a) = reshape(c1a,(/sys%Nocc_a*sys%Nunocc_a/))
                        qout(sys%Nocc_a*sys%Nunocc_a+1:sys%Nocc_a*sys%Nunocc_a+sys%Nocc_b*sys%Nunocc_b) = &
                        reshape(c1b,(/sys%Nocc_b*sys%Nunocc_b/))


                end subroutine update_C


end module davidson_module
