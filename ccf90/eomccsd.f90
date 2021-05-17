module eomccsd

        ! PRO TIP: If the EOMCC solver seems to "collapse" resulting in eigenvalues
        ! of 0, typically followed by a subsequent memory error, the culprit is 
        ! typically the initial guess. There could be 0's in places where there
        ! should not be, leading to trivial eigenvalue solution of the interaction
        ! matrix.

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use permutils, only: reorder_stripe
        use einsum_module, only: einsum
        use dgemm_module, only: gemm, gemmt, gemv, dot

        implicit none

        contains

               subroutine initial_guess(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,N,nroot,guess_type,B0,io)

                       use cis, only: cis_mat
                       use sort_module, only: argsort

                       integer, intent(in) :: io
                       type(sys_t), intent(in) :: sys
                       type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                       type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                       integer, intent(in) :: nroot, guess_type, N
                       real, intent(out) :: B0(N,nroot)
                       real, allocatable :: D1(:), C1(:,:), D2(:), w_cis(:), c1a(:), c1b(:)
                       integer, allocatable :: idx(:), idx2(:)
                       integer :: n1a, n1b, n2a, n2b, n2c, cnt, i, j, a, b

                       n1a = sys%Nunocc_a * sys%Nocc_a
                       n1b = sys%Nunocc_b * sys%Nocc_b
                       n2a = sys%Nunocc_a**2 * sys%Nocc_a**2
                       n2b = sys%Nunocc_a*sys%Nunocc_b*sys%Nocc_a*sys%Nocc_b
                       n2c = sys%Nunocc_b**2 * sys%Nocc_b**2

                       ! set initial guess to 0
                       B0 = 0.0

                       if (guess_type == 1) then ! DIAGONAL (KOOPMAN'S) GUESS

                          allocate(D1(n1a+n1b),idx(n1a+n1b))

                          cnt = 1
                          do i = 1,sys%Nocc_a
                             do a = 1,sys%Nunocc_a
                                D1(cnt) = H1A%uu(a,a)-H1A%oo(i,i)+H2A%uoou(a,i,i,a)
                                cnt = cnt + 1
                             end do
                          end do
                          do i = 1,sys%Nocc_b
                             do a = 1,sys%Nunocc_b
                                D1(cnt) = H1B%uu(a,a)-H1B%oo(i,i)+H2C%uoou(a,i,i,a)
                                cnt = cnt + 1
                             end do
                          end do

                        idx = argsort(D1)
                        do i = 1,nroot
                           B0(idx(i),i) = 1.0
                        end do

                        deallocate(D1,idx)

                     end if

                     if (guess_type == 2) then ! CIS GUESS

                       allocate(C1(n1a+n1b,n1a+n1b),w_cis(n1a+n1b),idx(n1a+n1b),c1a(n1a),c1b(n1b))

                       call cis_mat(sys,fA,fB,vA,vB,vC,C1,w_cis)
                       idx = argsort(w_cis)
                       w_cis = w_cis(idx)
                       C1 = C1(:,idx)

                       do i = 1,nroot
                          write(io,fmt=*) 'Root',i,' = ',w_cis(i),'HARTREE'
                          B0(1:n1a+n1b,i) = C1(:,i)
                       end do


                       deallocate(C1,w_cis,idx,c1a,c1b)

                    end if

                    if (guess_type == 3) then ! ONE-BODY + TWO-BODY DIAGONAL (KOOPMAN'S) GUESS

                          allocate(D1(n1a+n1b),idx(n1a+n1b))

                          cnt = 1
                          do i = 1,sys%Nocc_a
                             do a = 1,sys%Nunocc_a
                                D1(cnt) = H1A%uu(a,a)-H1A%oo(i,i)+H2A%uoou(a,i,i,a)
                                cnt = cnt + 1
                             end do
                          end do

                          do i = 1,sys%Nocc_b
                             do a = 1,sys%Nunocc_b
                                D1(cnt) = H1B%uu(a,a)-H1B%oo(i,i)+H2C%uoou(a,i,i,a)
                                cnt = cnt + 1
                             end do
                          end do

                        idx = argsort(D1)
                        do i = 1,nroot
                           B0(idx(i),i) = 1.0
                        end do
                        deallocate(D1,idx)

                        allocate(D2(n2a+n2b+n2c),idx2(n2a+n2b+n2c))

                        cnt = 1
                          do i = 1,sys%Nocc_a
                             do j = 1,sys%Nocc_a
                                do a = 1,sys%Nunocc_a
                                   do b = 1,sys%Nunocc_a
                                      D2(cnt) = H1A%uu(a,a)+H1A%uu(b,b)-H1A%oo(i,i)-H1A%oo(j,j)&
                                                +H2A%uoou(a,i,i,a)+H2A%uoou(a,j,j,a)+H2A%uoou(b,i,i,b)&
                                                +H2A%uoou(b,j,j,b)+H2A%uuuu(a,b,a,b)+H2A%oooo(i,j,i,j)
                                      cnt = cnt + 1
                                   end do
                                end do
                             end do
                          end do

                          do i = 1,sys%Nocc_a
                             do j = 1,sys%Nocc_b
                                do a = 1,sys%Nunocc_a
                                   do b = 1,sys%Nunocc_b
                                      D2(cnt) = H1A%uu(a,a)+H1B%uu(b,b)-H1A%oo(i,i)-H1B%oo(j,j)&
                                               +H2A%uoou(a,i,i,a)+H2C%uoou(b,j,j,b)-H2B%uouo(a,j,a,j)&
                                               -H2B%ouou(i,b,i,b)+H2B%uuuu(a,b,a,b)+H2B%oooo(i,j,i,j)
                                      cnt = cnt + 1
                                   end do
                                end do
                             end do
                          end do

                          do i = 1,sys%Nocc_b
                             do j = 1,sys%Nocc_b
                                do a = 1,sys%Nunocc_b
                                   do b = 1,sys%Nunocc_b
                                      D2(cnt) = H1B%uu(a,a)+H1B%uu(b,b)-H1B%oo(i,i)-H1B%oo(j,j)&
                                               +H2C%uoou(a,i,i,a)+H2C%uoou(a,j,j,a)+H2C%uoou(b,i,i,b)&
                                               +H2C%uoou(b,j,j,b)+H2C%oooo(i,j,i,j)+H2C%uuuu(a,b,a,b)
                                      cnt = cnt + 1
                                   end do
                                end do
                             end do
                          end do

                        idx2 = argsort(D2)
                        do i = 1,nroot
                           B0(n1a+n1b+idx2(i),i) = 1.0
                        end do

                        deallocate(D2,idx2)

                        call gramschmidt(B0)

                    end if

                    !if (initial_guess == 3.0) then ! CISd guess

                    !end if
                       

               end subroutine initial_guess
                       

               subroutine davidson_eomccsd(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                               t1a,t1b,t2a,t2b,t2c,N,nroot,VR,wR,tol,maxit,init_guess,io)

                        use dgeev_module, only: eig
                        use sort_module, only: argsort

                        integer, intent(in) :: io
                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        integer, intent(in) :: N
                        real, intent(in) :: tol
                        integer, intent(in) :: nroot, maxit, init_guess
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: VR(N,nroot), wR(nroot)
                        integer :: crsz, i, j, a, k, num_add, maxsz, it, cnt, idx(N), crsz0
                        real, allocatable :: sigma(:,:), aR(:,:), evR(:), evI(:), G(:,:),&
                                             rv(:), v(:), q(:), B(:,:), addB(:,:), overlaps(:) 
                        integer, allocatable :: id(:), id2(:)
                        real :: resid, resv(nroot), B0(N,nroot)
                        integer, parameter :: nvec = 10
                        real, parameter :: threshold = 1.0e-04, hatoev = 27.2114

                        ! print header
                        write(io,fmt=*) ''
                        write(io,fmt=*) '+++++++++++++++EOMCCSD ROUTINE - BLOCK DAVIDSON SOLVER++++++++++++++++'
                        write(io,fmt=*) ''
                        write(io,fmt=*) 'NUMBER OF ROOTS - ',nroot
                        write(io,fmt=*) 'MAX SUBSPACE SZ - ',nroot*nvec
                        write(io,fmt=*) 'CONVRG THRESH = ',tol
                        write(io,fmt=*) 'SUBSPACE EXPAND THRESH = ',threshold
                        if (init_guess == 1) then
                                write(io,fmt=*) 'INITIAL GUESS TYPE - DIAGONAL'
                        end if
                        if (init_guess == 2) then
                                write(io,fmt=*) 'INITIAL GUESS TYPE - CIS'
                        end if

                        ! Get the initial guess
                        call initial_guess(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,N,nroot,init_guess,B0,io)

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

                        ! find the projected subspace eigenvectors that correspond
                        ! to the initial guess
                        VR = B0

                        ! For the sake of the first iteration:
                        ! (1) initialize residuals to be large upon entering the Davidson loop
                        resv(1:nroot) = 100.0
                        ! (2) set num_add equal to nroot
                        num_add = nroot

                        ! main davidson loop
                        do it = 1,maxit

                                write(io,fmt=*) ''
                                write(io,'(1x,a13,1x,i0,10x,a15,1x,i0)') 'Iteration - ',it,'Subspace dim - ',crsz
                                write(io,fmt=*) '-----------------------------------------------------------'

                                ! allocate arrays
                                allocate(aR(crsz,crsz),evR(crsz),evI(crsZ),G(crsz,crsz),id(crsz),overlaps(crsz),id2(crsz))

                                ! check if subspace is stangated
                                if (crsz0 == crsz) then
                                   write(io,fmt=*) 'Davidson stagnated!'
                                   exit
                                else
                                   crsz0 = crsz
                                end if

                                ! orthogonalize guess space
                                call gramschmidt(B(:,1:crsz))

                                ! matrix-vector product
                                k = max(crsz-num_add+1,1)
                                call HR_matmat(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,&
                                        B(:,k:crsz),N,crsz-k+1,sigma(:,k:crsz))

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
                                      call update_R(rv,sys,H1A,H1B,evR(j),q)
                                      ! normalize update vector
                                      q = q/norm2(q)
                                      ! orthogonalize against all vectors in current search space
                                      call orth(q,B(:,1:crsz))
                                      ! and those that are to be added to the search space
                                      call orth(q,addB(:,1:num_add))

                                      ! add q to search space if its norm is large enough
                                      if (norm2(q) > threshold) then
                                          num_add = num_add + 1
                                          addB(:,num_add) = q
                                      end if

                                   end if

                                end do

                                ! orthonormalize the R vectors - important to enforce
                                ! biorthonormality in the left EOM solver
                                call gramschmidt(VR)

                                ! print status of each root
                                do j = 1,nroot
                                   write(io,'(1x,a8,i0,8x,a4,f12.8,8x,a6,f12.8)') 'Root - ',j,'e = ',wR(j),'|r| = ',resv(j)
                                end do

                                ! check exit condition
                                if (resid <= tol) then
                                   write(io,fmt='(/a)') 'EOMCCSD solver converged!'
                                   do i = 1,nroot
                                      write(io,'(1x,a4,2x,i0,a12,f12.8,2x,a10,a12,f12.8,2x,a2)') 'Root',i,'E = ',wR(i),'HARTREE',&
                                              'VEE =',wR(i)*hatoev,'eV'
                                   end do
                                   exit

                                ! otherwise, check if number of search vectors exceeds max size
                                elseif (crsz+num_add > maxsz) then

                                   ! collapse subspace    
                                   call gemm(B(:,1:crsz),aR(1:crsz,1:nroot),addB)
                                   B(:,1:nroot) = addB
                                   ! update search subspace dimension
                                   crsz = nroot

                                ! if not, expand the subspace
                                else

                                    ! add vectors to subspace
                                    B(:,crsz+1:crsz+num_add) = addB
                                    ! update search subspace dimension
                                    crsz = crsz + num_add

                                end if

                                ! deallocate temporary arrays
                                deallocate(aR,evR,evI,G,id,overlaps,id2)

                        end do

                        ! deallocate rest of the arrays
                        deallocate(rv,q,B,sigma,addB,v)

               end subroutine davidson_eomccsd

               subroutine update_R(q,sys,H1A,H1B,omega,qout)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: H1A, H1B
                        real, intent(in) :: q(:), omega
                        real, intent(out) :: qout(sys%Nocc_a*sys%Nunocc_a+sys%Nocc_b*sys%Nunocc_b + &
                                                  sys%Nunocc_a**2*sys%Nocc_a**2+sys%Nunocc_b**2*sys%Nocc_b**2 + &
                                                  sys%Nunocc_a*sys%Nunocc_b*sys%Nocc_a*sys%Nocc_b)
                        real :: r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                denom
                        integer :: a, b, i, j, n1a, n1b, n2a, n2b, n2c

                        n1a = sys%Nocc_a * sys%Nunocc_a
                        n1b = sys%Nocc_b * sys%Nunocc_b
                        n2a = sys%Nocc_a**2 * sys%Nunocc_a**2
                        n2b = sys%Nocc_a*sys%Nocc_b*sys%Nunocc_a*sys%Nunocc_b
                        n2c = sys%Nocc_b**2 * sys%Nunocc_b**2

                        r1a = reshape(q(1:n1a),shape(r1a))
                        r1b = reshape(q(n1a+1:n1a+n1b),shape(r1b))
                        r2a = reshape(q(n1b+n1a+1:n1a+n1b+n2a),shape(r2a))
                        r2b = reshape(q(n2a+n1b+n1a+1:n1a+n1b+n2a+n2b),shape(r2b))
                        r2c = reshape(q(n2b+n2a+n1b+n1a+1:n1a+n1b+n2a+n2b+n2c),shape(r2c))

                        do i = 1,sys%Nocc_a
                           do a = 1,sys%Nunocc_a
                              denom = H1A%uu(a,a)-H1A%oo(i,i)
                              r1a(a,i) = r1a(a,i)/(omega-denom)
                           end do
                        end do

                        do i = 1,sys%Nocc_b
                           do a = 1,sys%Nunocc_b
                              denom = H1B%uu(a,a)-H1B%oo(i,i)
                              r1b(a,i) = r1b(a,i)/(omega-denom)
                           end do
                        end do

                        do i = 1,sys%Nocc_a
                           do j = 1,sys%Nocc_a
                              do a = 1,sys%Nunocc_a
                                 do b = 1,sys%Nunocc_a
                                    denom = H1A%uu(a,a)+H1A%uu(b,b)-H1A%oo(i,i)-H1A%oo(j,j)
                                    r2a(b,a,j,i) = r2a(b,a,j,i)/(omega-denom)
                                    !r2a(a,b,j,i) = -r2a(b,a,j,i)
                                    !r2a(b,a,i,j) = -r2a(b,a,j,i)
                                    !r2a(a,b,i,j) = r2a(b,a,j,i)
                                 end do
                              end do
                           end do
                        end do

                        do j = 1,sys%Nocc_b
                           do i = 1,sys%Nocc_a
                              do b = 1,sys%Nunocc_b
                                 do a = 1,sys%Nunocc_a
                                    denom = H1A%uu(a,a)+H1B%uu(b,b)-H1A%oo(i,i)-H1B%oo(j,j)
                                    r2b(a,b,i,j) = r2b(a,b,i,j)/(omega-denom)
                                 end do
                              end do
                           end do
                        end do

                        do i = 1,sys%Nocc_b
                           do j = 1,sys%Nocc_b
                              do a = 1,sys%Nunocc_b
                                 do b = 1,sys%Nunocc_b
                                    denom = H1B%uu(a,a)+H1B%uu(b,b)-H1B%oo(i,i)-H1B%oo(j,j)
                                    r2c(b,a,j,i) = r2c(b,a,j,i)/(omega-denom)
                                    !r2c(a,b,j,i) = -r2c(b,a,j,i)
                                    !r2c(b,a,i,j) = -r2c(b,a,j,i)
                                    !r2c(a,b,i,j) = r2c(b,a,j,i)
                                 end do
                              end do
                           end do
                        end do

                        qout(1:n1a) = reshape(r1a,(/n1a/))
                        qout(n1a+1:n1a+n1b) = reshape(r1b,(/n1b/))
                        qout(n1b+n1a+1:n1a+n1b+n2a) = reshape(r2a,(/n2a/))
                        qout(n2a+n1b+n1a+1:n1a+n1b+n2a+n2b) = reshape(r2b,(/n2b/))
                        qout(n2b+n2a+n1b+n1a+1:n1a+n1b+n2a+n2b+n2c) = reshape(r2c,(/n2c/))


               end subroutine update_R

               subroutine HR_matvec(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,Rvec,ndim,sigma)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        integer, intent(in) :: ndim
                        real, intent(in) :: Rvec(ndim)
                        real, intent(out) :: sigma(ndim)
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real :: r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                X1A(sys%Nunocc_a,sys%Nocc_a), X1B(sys%Nunocc_b,sys%Nocc_b), &
                                X2A(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                X2B(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                X2C(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                sigma_1a(sys%Nunocc_a*sys%Nocc_a), sigma_1b(sys%Nunocc_b*sys%Nocc_b), &
                                sigma_2a(sys%Nunocc_a**2*sys%Nocc_a**2), &
                                sigma_2b(sys%Nunocc_a*sys%Nunocc_b*sys%Nocc_a*sys%Nocc_b), &
                                sigma_2c(sys%Nunocc_b**2*sys%Nocc_b**2)
                        integer :: n1a, n1b, n2a, n2b, n2c

                        n1a = sys%Nunocc_a * sys%Nocc_a
                        n1b = sys%Nunocc_b * sys%Nocc_b
                        n2a = sys%Nunocc_a**2 * sys%Nocc_a**2
                        n2b = sys%Nunocc_a * sys%Nunocc_b * sys%Nocc_a * sys%Nocc_b
                        n2c = sys%Nunocc_b**2 * sys%Nocc_b**2

                        r1a = reshape(Rvec(1:n1a),shape(r1a))
                        r1b = reshape(Rvec(n1a+1:n1a+n1b),shape(r1b))
                        r2a = reshape(Rvec(n1b+n1a+1:n1a+n1b+n2a),shape(r2a))
                        r2b = reshape(Rvec(n2a+n1b+n1a+1:n1a+n1b+n2a+n2b),shape(r2b))
                        r2c = reshape(Rvec(n2b+n2a+n1b+n1a+1:n1a+n1b+n2a+n2b+n2c),shape(r2c))

                        call HR_1A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                r1a,r1b,r2a,r2b,r2c,X1A)
                        sigma_1a = reshape(X1A,(/n1a/))
                        call HR_1B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                r1a,r1b,r2a,r2b,r2c,X1B)
                        sigma_1b = reshape(X1B,(/n1b/))
                        call HR_2A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                r1a,r1b,r2a,r2b,r2c,X2A)
                        sigma_2a = reshape(X2A,(/n2a/))
                        call HR_2B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                r1a,r1b,r2a,r2b,r2c,X2B)
                        sigma_2b = reshape(X2B,(/n2b/))
                        call HR_2C(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                r1a,r1b,r2a,r2b,r2c,X2C)
                        sigma_2c = reshape(X2C,(/n2c/))

                        sigma(1:n1a) = sigma_1a
                        sigma(n1a+1:n1a+n1b) = sigma_1b
                        sigma(n1a+n1b+1:n1a+n1b+n2a) = sigma_2a
                        sigma(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b) = sigma_2b
                        sigma(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c) = sigma_2c


                end subroutine HR_matvec

               subroutine HR_matmat(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,Rvec,ndim,crsz,sigma)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        integer, intent(in) :: ndim, crsz
                        real, intent(in) :: Rvec(ndim,crsz)
                        real, intent(out) :: sigma(ndim,crsz)
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real :: r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                X1A(sys%Nunocc_a,sys%Nocc_a), X1B(sys%Nunocc_b,sys%Nocc_b), &
                                X2A(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                X2B(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                X2C(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                sigma_1a(sys%Nunocc_a*sys%Nocc_a), sigma_1b(sys%Nunocc_b*sys%Nocc_b), &
                                sigma_2a(sys%Nunocc_a**2*sys%Nocc_a**2), &
                                sigma_2b(sys%Nunocc_a*sys%Nunocc_b*sys%Nocc_a*sys%Nocc_b), &
                                sigma_2c(sys%Nunocc_b**2*sys%Nocc_b**2)
                        integer :: j, n1a, n1b, n2a, n2b, n2c

                        n1a = sys%Nunocc_a * sys%Nocc_a
                        n1b = sys%Nunocc_b * sys%Nocc_b
                        n2a = sys%Nunocc_a**2 * sys%Nocc_a**2
                        n2b = sys%Nunocc_a * sys%Nunocc_b * sys%Nocc_a * sys%Nocc_b
                        n2c = sys%Nunocc_b**2 * sys%Nocc_b**2

                        do j = 1,crsz

                                r1a = reshape(Rvec(1:n1a,j),shape(r1a))
                                r1b = reshape(Rvec(n1a+1:n1a+n1b,j),shape(r1b))
                                r2a = reshape(Rvec(n1b+n1a+1:n1a+n1b+n2a,j),shape(r2a))
                                r2b = reshape(Rvec(n2a+n1b+n1a+1:n1a+n1b+n2a+n2b,j),shape(r2b))
                                r2c = reshape(Rvec(n2b+n2a+n1b+n1a+1:n1a+n1b+n2a+n2b+n2c,j),shape(r2c))

                                call HR_1A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                        r1a,r1b,r2a,r2b,r2c,X1A)
                                sigma_1a = reshape(X1A,(/n1a/))
                                call HR_1B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                        r1a,r1b,r2a,r2b,r2c,X1B)
                                sigma_1b = reshape(X1B,(/n1b/))
                                call HR_2A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                        r1a,r1b,r2a,r2b,r2c,X2A)
                                sigma_2a = reshape(X2A,(/n2a/))
                                call HR_2B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                        r1a,r1b,r2a,r2b,r2c,X2B)
                                sigma_2b = reshape(X2B,(/n2b/))
                                call HR_2C(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c, &
                                        r1a,r1b,r2a,r2b,r2c,X2C)
                                sigma_2c = reshape(X2C,(/n2c/))

                                sigma(1:n1a,j) = sigma_1a
                                sigma(n1a+1:n1a+n1b,j) = sigma_1b
                                sigma(n1a+n1b+1:n1a+n1b+n2a,j) = sigma_2a
                                sigma(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b,j) = sigma_2b
                                sigma(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c,j) = sigma_2c

                        end do

                end subroutine HR_matmat

                subroutine HR_1A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                            r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: HR(sys%Nunocc_a,sys%Nocc_a)
                        real, allocatable :: Z(:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                allocate(Z(nua,noa))
                                call einsum('mi,am->ai',H1A%oo,r1a,Z)
                                HR = -1.0*Z
                                call einsum('ae,ei->ai',H1A%uu,r1a,Z)
                                HR = HR + Z
                                call einsum('amie,em->ai',H2A%uoou,r1a,Z)
                                HR = HR + Z
                                call einsum('amie,em->ai',H2B%uoou,r1b,Z)
                                HR = HR + Z
                                call einsum('mnif,afmn->ai',H2A%ooou,r2a,Z)
                                HR = HR - 0.5*Z
                                call einsum('mnif,afmn->ai',H2B%ooou,r2b,Z)
                                HR = HR - Z
                                call einsum('anef,efin-.ai',H2A%uouu,r2a,Z)
                                HR = HR + 0.5*Z
                                call einsum('anef,efin->ai',H2B%uouu,r2b,Z)
                                HR = HR + Z
                                call einsum('me,aeim->ai',H1A%ou,r2a,Z)
                                HR = HR + Z
                                call einsum('me,aeim->ai',H1B%ou,r2b,Z)
                                HR = HR + Z
                                deallocate(Z)

                        end associate

                end subroutine HR_1A

                subroutine HR_1B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                            r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: HR(sys%Nunocc_b,sys%Nocc_b)
                        real, allocatable :: Z(:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                allocate(Z(nub,nob))
                                call einsum('mi,am->ai',H1B%oo,r1b,Z)
                                HR = -1.0*Z
                                call einsum('ae,ei->ai',H1B%uu,r1b,Z)
                                HR = HR + Z
                                call einsum('amie,em->ai',H2C%uoou,r1b,Z)
                                HR = HR + Z
                                call einsum('maei,em->ai',H2B%ouuo,r1a,Z)
                                HR = HR + Z
                                call einsum('mnif,afmn->ai',H2C%ooou,r2c,Z)
                                HR = HR - 0.5*Z
                                call einsum('nmfi,fanm->ai',H2B%oouo,r2b,Z)
                                HR = HR - Z
                                call einsum('anef,efin-.ai',H2C%uouu,r2c,Z)
                                HR = HR + 0.5*Z
                                call einsum('nafe,feni->ai',H2B%ouuu,r2b,Z)
                                HR = HR + Z
                                call einsum('me,aeim->ai',H1B%ou,r2c,Z)
                                HR = HR + Z
                                call einsum('me,eami->ai',H1A%ou,r2b,Z)
                                HR = HR + Z
                                deallocate(Z)

                        end associate


                end subroutine HR_1B
        

                subroutine HR_2A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                            r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: HR(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a)
                        real, allocatable :: Z(:,:), Q(:,:,:,:), D1(:,:,:,:), D2(:,:,:,:), D3(:,:,:,:), &
                                             D4(:,:,:,:), D5(:,:,:,:), D6(:,:,:,:), D7(:,:,:,:), D8(:,:,:,:), &
                                             D9(:,:,:,:), D10(:,:,:,:), D11(:,:,:,:), D12(:,:,:,:), D13(:,:,:,:), &
                                             D14(:,:,:,:), D15(:,:,:,:), D16(:,:,:,:), &
                                             Dij(:,:,:,:), Dab(:,:,:,:), Dabij(:,:,:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                allocate(D1(nua,nua,noa,noa),D2(nua,nua,noa,noa),D3(nua,nua,noa,noa),D4(nua,nua,noa,noa),&
                                         D5(nua,nua,noa,noa),D6(nua,nua,noa,noa),D7(nua,nua,noa,noa),D8(nua,nua,noa,noa),&
                                         D9(nua,nua,noa,noa),D10(nua,nua,noa,noa),D11(nua,nua,noa,noa),D12(nua,nua,noa,noa),&
                                         D13(nua,nua,noa,noa),D14(nua,nua,noa,noa),D15(nua,nua,noa,noa),D16(nua,nua,noa,noa))

                                call einsum('mi,abmj->abij',-1.0*H1A%oo,r2a,D1)
                                call einsum('ae,ebij->abij',H1A%uu,r2a,D2)
                                
                                allocate(Q(nua,nua,noa,noa))
                                call einsum('mnij,abmn->abij',H2A%oooo,r2a,Q)
                                HR = 0.5*Q
                                call einsum('abef,efij->abij',H2A%uuuu,r2a,Q)
                                HR = HR + 0.5*Q
                                deallocate(Q)

                                call einsum('amie,ebmj->abij',H2A%uoou,r2a,D3)
                                call einsum('amie,bejm->abij',H2B%uoou,r2b,D4)
                                call einsum('bmji,am->abij',-1.0*H2A%uooo,r1a,D5)
                                call einsum('baje,ei->abij',H2A%uuou,r1a,D6)

                                allocate(Z(nua,nua))
                                call einsum('mnef,bfmn->eb',-0.5*vA%oouu,r2a,Z)
                                call einsum('eb,aeij->abij',Z,t2a,D7)
                                call einsum('mnef,bfmn->eb',-1.0*vB%oouu,r2b,Z)
                                call einsum('eb,aeij->abij',Z,t2a,D8)
                                deallocate(Z)

                                allocate(Z(noa,noa))
                                call einsum('mnef,efjn->mj',0.5*vA%oouu,r2a,Z)
                                call einsum('mj,abim->abij',-1.0*Z,t2a,D9)
                                call einsum('mnef,efjn->mj',vB%oouu,r2b,Z)
                                call einsum('mj,abim->abij',-1.0*Z,t2a,D10)
                                deallocate(Z)

                                allocate(Z(nua,nua))
                                call einsum('amfe,em->af',H2A%uouu,r1a,Z)
                                call einsum('af,fbij->abij',Z,t2a,D11)
                                deallocate(Z)

                                allocate(Z(noa,noa))
                                call einsum('nmie,em->ni',H2A%ooou,r1a,Z)
                                call einsum('ni,abnj->abij',-1.0*Z,t2a,D12)
                                deallocate(Z)

                                allocate(Z(nua,nua))
                                call einsum('amfe,em->af',H2B%uouu,r1b,Z)
                                call einsum('af,fbij->abij',Z,t2a,D13)
                                deallocate(Z)

                                allocate(Z(noa,noa))
                                call einsum('nmie,em->ni',H2B%ooou,r1b,Z)
                                call einsum('ni,abnj->abij',-1.0*Z,t2a,D14)
                                deallocate(Z)

                                allocate(Dij(nua,nua,noa,noa),Dab(nua,nua,noa,noa),Dabij(nua,nua,noa,noa))

                                Dij = D1 + D6 + D9 + D10 + D12 + D14
                                Dab = D2 + D5 + D7 + D8 + D11 + D13
                                Dabij = D3 + D4

                                call reorder_stripe(4,shape(Dij),size(Dij),'1243',Dij,D1)
                                Dij = Dij - D1
                                call reorder_stripe(4,shape(Dab),size(Dab),'2134',Dab,D2)
                                Dab = Dab - D2
                                call reorder_stripe(4,shape(Dabij),size(Dabij),'2134',Dabij,D1)
                                call reorder_stripe(4,shape(Dabij),size(Dabij),'1243',Dabij,D2)
                                call reorder_stripe(4,shape(Dabij),size(Dabij),'2143',Dabij,D3)
                                Dabij = Dabij - D1 - D2 + D3

                                HR = HR + Dij + Dab + Dabij

                                deallocate(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,Dij,Dab,Dabij)
                                
                        end associate

                end subroutine HR_2A

                subroutine HR_2B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                            r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: HR(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b)
                        real, allocatable :: Z(:,:), Z1(:,:), Q(:,:,:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                allocate(Q(nua,nub,noa,nob))
                                call einsum('ae,ebij->abij',H1A%uu,r2b,Q)
                                HR = Q
                                call einsum('be,aeij->abij',H1B%uu,r2b,Q)
                                HR = HR + Q
                                call einsum('mi,abmj->abij',H1A%oo,r2b,Q)
                                HR = HR - Q
                                call einsum('mj,abim->abij',H1B%oo,r2b,Q)
                                HR = HR - Q
                                call einsum('mnij,abmn->abij',H2B%oooo,r2b,Q)
                                HR = HR + Q
                                call einsum('abef,efij->abij',H2B%uuuu,r2b,Q)
                                HR = HR + Q
                                call einsum('amie,ebmj->abij',H2A%uoou,r2b,Q)
                                HR = HR + Q
                                call einsum('amie,ebmj->abij',H2B%uoou,r2c,Q)
                                HR = HR + Q
                                call einsum('mbej,aeim->abij',H2B%ouuo,r2a,Q)
                                HR = HR + Q
                                call einsum('bmje,aeim->abij',H2C%uoou,r2b,Q)
                                HR = HR + Q
                                call einsum('mbie,aemj->abij',H2B%ouou,r2b,Q)
                                HR = HR - Q
                                call einsum('amej,ebim->abij',H2B%uouo,r2b,Q)
                                HR = HR - Q
                                call einsum('abej,ei->abij',H2B%uuuo,r1a,Q)
                                HR = HR + Q
                                call einsum('abie,ej->abij',H2B%uuou,r1b,Q)
                                HR = HR + Q
                                call einsum('mbij,am->abij',H2B%ouoo,r1a,Q)
                                HR = HR - Q
                                call einsum('amij,bm->abij',H2B%uooo,r1b,Q)
                                HR = HR - Q
                                
                                allocate(Z(nua,nua),Z1(nua,nua))
                                call einsum('mnef,afmn->ae',-0.5*vA%oouu,r2a,Z)
                                call einsum('mnef,afmn->ae',-1.0*vB%oouu,r2b,Z1)
                                Z = Z + Z1
                                call einsum('ae,ebij->abij',Z,t2b,Q)
                                HR = HR + Q
                                call einsum('amfe,em->af',H2A%uouu,r1a,Z)
                                call einsum('amfe,em->af',H2B%uouu,r1b,Z1)
                                Z = Z + Z1
                                call einsum('af,fbij->abij',Z,t2b,Q)
                                HR = HR + Q
                                deallocate(Z,Z1)

                                allocate(Z(noa,noa),Z1(noa,noa))
                                call einsum('mnef,efin->mi',0.5*vA%oouu,r2a,Z)
                                call einsum('mnef,efin->mi',vB%oouu,r2b,Z1)
                                Z = Z + Z1
                                call einsum('mi,abmj->abij',Z,t2b,Q)
                                HR = HR - Q
                                call einsum('nmie,em->ni',H2A%ooou,r1a,Z)
                                call einsum('nmie,em->ni',H2B%ooou,r1b,Z1)
                                Z = Z + Z1
                                call einsum('ni,abnj->abij',Z,t2b,Q)
                                HR = HR - Q
                                deallocate(Z,Z1)

                                allocate(Z(nub,nub),Z1(nub,nub))
                                call einsum('nmfe,fbnm->be',-1.0*vB%oouu,r2b,Z)
                                call einsum('mnef,bfmn->be',-0.5*vC%oouu,r2c,Z1)
                                Z = Z + Z1
                                call einsum('be,aeij->abij',Z,t2b,Q)
                                HR = HR + Q
                                call einsum('mbef,em->bf',H2B%ouuu,r1a,Z)
                                call einsum('bmfe,em->bf',H2C%uouu,r1b,Z1)
                                Z = Z + Z1
                                call einsum('bf,afij->abij',Z,t2b,Q)
                                HR = HR + Q                              
                                deallocate(Z,Z1)

                                allocate(Z(nob,nob),Z1(nob,nob))
                                call einsum('nmfe,fenj->mj',vB%oouu,r2b,Z)
                                call einsum('mnef,efjn->mj',0.5*vC%oouu,r2c,Z1)
                                Z = Z + Z1
                                call einsum('mj,abim->abij',Z,t2b,Q)
                                HR = HR - Q
                                call einsum('mnej,em->nj',H2B%oouo,r1a,Z)
                                call einsum('nmje,em->nj',H2C%ooou,r1b,Z1)
                                Z = Z + Z1
                                call einsum('nj,abin->abij',Z,t2b,Q)
                                HR = HR - Q
                                deallocate(Q,Z,Z1)


                        end associate


                end subroutine HR_2B

                subroutine HR_2C(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                            r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: HR(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, allocatable :: Z(:,:), Q(:,:,:,:), D1(:,:,:,:), D2(:,:,:,:), D3(:,:,:,:), &
                                             D4(:,:,:,:), D5(:,:,:,:), D6(:,:,:,:), D7(:,:,:,:), D8(:,:,:,:), &
                                             D9(:,:,:,:), D10(:,:,:,:), D11(:,:,:,:), D12(:,:,:,:), D13(:,:,:,:), &
                                             D14(:,:,:,:), D15(:,:,:,:), D16(:,:,:,:), &
                                             Dij(:,:,:,:), Dab(:,:,:,:), Dabij(:,:,:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                allocate(D1(nub,nub,nob,nob),D2(nub,nub,nob,nob),D3(nub,nub,nob,nob),D4(nub,nub,nob,nob),&
                                         D5(nub,nub,nob,nob),D6(nub,nub,nob,nob),D7(nub,nub,nob,nob),D8(nub,nub,nob,nob),&
                                         D9(nub,nub,nob,nob),D10(nub,nub,nob,nob),D11(nub,nub,nob,nob),D12(nub,nub,nob,nob),&
                                         D13(nub,nub,nob,nob),D14(nub,nub,nob,nob),D15(nub,nub,nob,nob),D16(nub,nub,nob,nob))

                                call einsum('mi,abmj->abij',-1.0*H1B%oo,r2c,D1)
                                call einsum('ae,ebij->abij',H1B%uu,r2c,D2)
                                
                                allocate(Q(nub,nub,nob,nob))
                                call einsum('mnij,abmn->abij',H2C%oooo,r2c,Q)
                                HR = 0.5*Q
                                call einsum('abef,efij->abij',H2C%uuuu,r2c,Q)
                                HR = HR + 0.5*Q
                                deallocate(Q)

                                call einsum('amie,ebmj->abij',H2C%uoou,r2c,D3)
                                call einsum('maei,ebmj->abij',H2B%ouuo,r2b,D4)
                                call einsum('bmji,am->abij',-1.0*H2C%uooo,r1b,D5)
                                call einsum('baje,ei->abij',H2C%uuou,r1b,D6)

                                allocate(Z(nub,nub))
                                call einsum('mnef,bfmn->eb',-0.5*vC%oouu,r2c,Z)
                                call einsum('eb,aeij->abij',Z,t2c,D7)
                                call einsum('nmfe,fbnm->eb',-1.0*vB%oouu,r2b,Z)
                                call einsum('eb,aeij->abij',Z,t2c,D8)
                                deallocate(Z)

                                allocate(Z(nob,nob))
                                call einsum('mnef,efjn->mj',0.5*vC%oouu,r2c,Z)
                                call einsum('mj,abim->abij',-1.0*Z,t2c,D9)
                                call einsum('nmfe,fenj->mj',vB%oouu,r2b,Z)
                                call einsum('mj,abim->abij',-1.0*Z,t2c,D10)
                                deallocate(Z)

                                allocate(Z(nub,nub))
                                call einsum('amfe,em->af',H2C%uouu,r1b,Z)
                                call einsum('af,fbij->abij',Z,t2c,D11)
                                deallocate(Z)

                                allocate(Z(nob,nob))
                                call einsum('nmie,em->ni',H2C%ooou,r1b,Z)
                                call einsum('ni,abnj->abij',-1.0*Z,t2c,D12)
                                deallocate(Z)

                                allocate(Z(nub,nub))
                                call einsum('maef,em->af',H2B%ouuu,r1a,Z)
                                call einsum('af,fbij->abij',Z,t2c,D13)
                                deallocate(Z)

                                allocate(Z(noa,noa))
                                call einsum('mnei,em->ni',H2B%oouo,r1a,Z)
                                call einsum('ni,abnj->abij',-1.0*Z,t2c,D14)
                                deallocate(Z)

                                allocate(Dij(nub,nub,nob,nob),Dab(nub,nub,nob,nob),Dabij(nub,nub,nob,nob))

                                Dij = D1 + D6 + D9 + D10 + D12 + D14
                                Dab = D2 + D5 + D7 + D8 + D11 + D13
                                Dabij = D3 + D4

                                call reorder_stripe(4,shape(Dij),size(Dij),'1243',Dij,D1)
                                Dij = Dij - D1
                                call reorder_stripe(4,shape(Dab),size(Dab),'2134',Dab,D2)
                                Dab = Dab - D2
                                call reorder_stripe(4,shape(Dabij),size(Dabij),'2134',Dabij,D1)
                                call reorder_stripe(4,shape(Dabij),size(Dabij),'1243',Dabij,D2)
                                call reorder_stripe(4,shape(Dabij),size(Dabij),'2143',Dabij,D3)
                                Dabij = Dabij - D1 - D2 + D3

                                HR = HR + Dij + Dab + Dabij

                                deallocate(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,Dij,Dab,Dabij)
                                
                        end associate


                end subroutine HR_2C

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

                subroutine orth(q,B)

                        real, intent(inout) :: q(:)
                        real, intent(in) :: B(:,:)
                        real, allocatable :: v(:)
                        integer :: i, n
                        real :: s

                        n = size(B,2)

                        allocate(v(n))
                        do i = 1,n
                           v = B(:,i)/norm2(B(:,i))
                           s = dot(v,q)
                           q = q - s*v
                        end do
                        deallocate(v)

                end subroutine orth

end module eomccsd
