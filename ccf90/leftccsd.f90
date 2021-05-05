module leftccsd

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use permutils, only: reorder_stripe
        use einsum_module, only: einsum
        use dgemm_module, only: gemm, gemmt, gemv, dot

        implicit none

        contains

                subroutine left_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,L,R,&
                                H1A,H1B,H2A,H2B,H2C,tol,ndiis,maxit,shift,iroot,omega)

                        use diis, only: diis_xtrap
                        use cc_energy, only: calc_cc_energy
                
                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            tol, shift, omega(:), R(:,:)
                        integer, intent(in) :: maxit, ndiis, iroot
                        real :: l1a(sys%Nunocc_a,sys%Nocc_a), l1b(sys%Nunocc_b,sys%Nocc_b), &
                                l2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                l2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                l2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: L(sys%Nocc_a**2*sys%Nunocc_a**2+sys%Nocc_b**2*sys%Nunocc_b**2+&
                                                  sys%Nocc_a*sys%Nocc_b*sys%Nunocc_a*sys%Nunocc_b+&
                                                  sys%Nocc_a*sys%Nunocc_a+sys%Nocc_b*sys%Nunocc_b)
                        integer :: n1a, n1b, n2a, n2b, n2c, ndim, it, it_diis, noa, nob, nua, nub, i, nroot
                        real, allocatable :: L_resid_list(:,:), L_list(:,:), X1A(:,:), X1B(:,:), &
                                             X2A(:,:,:,:), X2B(:,:,:,:), X2C(:,:,:,:), LH(:), L_resid(:)
                        real :: resid, Ecorr, Elcc, e0, LR, val

                        write(*,fmt=*) ''
                        write(*,fmt=*) '++++++++++++++++++++ LEFT CCSD ROUTINE +++++++++++++++++++++'
                        write(*,fmt=*) ''

                        noa = sys%Nocc_a
                        nob = sys%Nocc_b
                        nua = sys%Nunocc_a
                        nub = sys%Nunocc_b

                        n1a = sys%Nocc_a * sys%Nunocc_a
                        n1b = sys%Nocc_b * sys%Nunocc_b
                        n2a = sys%Nunocc_a**2 * sys%Nocc_a**2
                        n2b = sys%Nunocc_a*sys%Nunocc_b*sys%Nocc_a*sys%Nocc_b
                        n2c = sys%Nunocc_b**2 * sys%Nocc_b**2
                        ndim = n1a + n1b + n2a + n2b + n2c

                        nroot = size(R,2)

                        call calc_cc_energy(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Ecorr)

                        if (iroot == 0) then
                                e0 = 0.0
                                L(1:n1a) = reshape(t1a,(/n1a/))
                                L(n1a+1:n1a+n1b) = reshape(t1b,(/n1b/))
                                L(n1a+n1b+1:n1a+n1b+n2a) = reshape(t2a,(/n2a/))
                                L(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b) = reshape(t2b,(/n2b/))
                                L(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c) = reshape(t2c,(/n2c/))
                        else
                                e0 = omega(iroot)
                                L = R(:,iroot)
                        end if

                        allocate(L_list(ndim,ndiis),L_resid_list(ndim,ndiis),X1A(nua,noa),X1B(nub,nob),&
                                X2A(nua,nua,noa,noa),X2B(nua,nub,noa,nob),X2C(nub,nub,nob,nob),LH(ndim),&
                                L_resid(ndim))

                        write(*,fmt=*) '               Iteration      E(corr)                       Residuum'         
                        write(*,fmt=*) '------------------------------------------------------------------------'

                        it_diis = 0

                        do it = 0,maxit

                           l1a = reshape(L(1:n1a),(/nua,noa/))
                           l1b = reshape(L(n1a+1:n1a+n1b),(/nub,nob/))
                           l2a = reshape(L(n1a+n1b+1:n1a+n1b+n2a),(/nua,nua,noa,noa/))
                           l2b = reshape(L(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b),(/nua,nub,noa,nob/))
                           l2c = reshape(L(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c),(/nub,nub,nob,nob/))
                
                           call LH_1A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,X1A)                           
                           call LH_1B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,X1B)
                           call LH_2A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,X2A)
                           call LH_2B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,X2B)
                           call LH_2C(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,X2C)

                           LH(1:n1a) = reshape(X1A,(/n1a/))
                           LH(n1a+1:n1a+n1b) = reshape(X1B,(/n1b/))
                           LH(n1a+n1b+1:n1a+n1b+n2a) = reshape(X2A,(/n2a/))
                           LH(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b) = reshape(X2B,(/n2b/))
                           LH(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c) = reshape(X2C,(/n2c/))

                           X1A = X1A - e0*l1a
                           X1B = X1B - e0*l1b
                           X2A = X2A - e0*l2a
                           X2B = X2B - e0*l2b
                           X2C = X2C - e0*l2c

                           call update_L(sys,l1a,l1b,l2a,l2b,l2c,H1A,H1B,X1A,X1B,X2A,X2B,X2C,shift,e0)

                           L(1:n1a) = reshape(l1a,(/n1a/))
                           L(n1a+1:n1a+n1b) = reshape(l1b,(/n1b/))
                           L(n1a+n1b+1:n1a+n1b+n2a) = reshape(l2a,(/n2a/))
                           L(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b) = reshape(l2b,(/n2b/))
                           L(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c) = reshape(l2c,(/n2c/))

                           L_resid = LH - e0*L
                           resid = dot(L_resid,L_resid)
                           resid = sqrt(resid)

                           if (iroot /= 0) then
                                LR = dot(L,R(:,iroot))
                                L = L/LR
                           end if

                           Elcc = norm2(LH)/norm2(L)

                           if (resid < tol) then

                                   write(*,fmt=*) 'Left CC converged!'

                                   if (iroot /= 0) then
                                        LR = 0.0
                                        do i = 1,nroot
                                           val = dot(L,R(:,i))
                                           LR = LR + val
                                        end do
                                        write(*,fmt=*) 'LR = ',LR
                                   end if

                                   exit
                           end if

                           L_list(:,mod(it,ndiis)+1) = L
                           L_resid_list(:,mod(it,ndiis)+1) = L_resid

                           if ( (mod(it,ndiis) == 0) .AND. it > 1) then
                              it_diis = it_diis + 1
                              write(*,fmt='(1x,a15,1x,i0)') 'DIIS cycle - ',it_diis
                              call diis_xtrap(L,L_list,L_resid_list)
                           end if
                                    
                           write(*,fmt=*) it,' ',Elcc,' ',resid

                       end do                       

                       deallocate(L_list,L_resid_list,X1A,X1B,X2A,X2B,X2C,LH,L_resid)
       
                end subroutine left_ccsd


                subroutine update_L(sys,l1a,l1b,l2a,l2b,l2c,H1A,H1B,X1A,X1B,X2A,X2B,X2C,shift,e0)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: H1A, H1B
                        real, intent(inout) :: l1a(sys%Nunocc_a,sys%Nocc_a), l1b(sys%Nunocc_b,sys%Nocc_b), &
                                               l2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                               l2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                               l2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(in) :: X1A(sys%Nunocc_a,sys%Nocc_a), X1B(sys%Nunocc_b,sys%Nocc_b), &
                                            X2A(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            X2B(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            X2C(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(in) :: shift, e0
                        real :: denom
                        integer :: i, j, a, b

                        do i = 1,sys%Nocc_a
                           do a = 1,sys%Nunocc_a
                              denom = H1A%uu(a,a) - H1A%oo(i,i)
                              l1a(a,i) = l1a(a,i) - X1A(a,i)/(denom+shift-e0)
                           end do
                        end do


                        do i = 1,sys%Nocc_b
                           do a = 1,sys%Nunocc_b
                            denom = H1B%uu(a,a) - H1B%oo(i,i)
                            l1b(a,i) = l1b(a,i) - X1B(a,i)/(denom+shift-e0)
                           end do
                        end do

                        do i = 1,sys%Nocc_a
                           do j = i+1,sys%Nocc_a
                              do a = 1,sys%Nunocc_a
                                 do b = a+1,sys%Nunocc_a
                                    denom = H1A%uu(a,a) + H1A%uu(b,b) - H1A%oo(i,i) - H1A%oo(j,j)
                                    l2a(b,a,j,i) = l2a(b,a,j,i) - X2A(b,a,j,i)/(denom+shift-e0)
                                    l2a(a,b,j,i) = -l2a(b,a,j,i)
                                    l2a(b,a,i,j) = -l2a(b,a,j,i)
                                    l2a(a,b,i,j) = l2a(b,a,j,i)
                                 end do
                              end do
                           end do
                        end do

                        do j = 1,sys%Nocc_b
                           do i = 1,sys%Nocc_a
                              do b = 1,sys%Nunocc_b
                                 do a = 1,sys%Nunocc_a
                                    denom = H1A%uu(a,a) + H1B%uu(b,b) - H1A%oo(i,i) - H1B%oo(j,j)
                                    l2b(a,b,i,j) = l2b(a,b,i,j) - X2B(a,b,i,j)/(denom+shift-e0)
                                 end do
                              end do
                           end do
                        end do

                        do i = 1,sys%Nocc_b
                           do j = i+1,sys%Nocc_b
                              do a = 1,sys%Nunocc_b
                                 do b = a+1,sys%Nunocc_b
                                    denom = H1B%uu(a,a) + H1B%uu(b,b) - H1B%oo(i,i) - H1B%oo(j,j)
                                    l2c(b,a,j,i) = l2c(b,a,j,i) - X2C(b,a,j,i)/(denom+shift-e0)
                                    l2c(a,b,j,i) = -l2c(b,a,j,i)
                                    l2c(b,a,i,j) = -l2c(b,a,j,i)
                                    l2c(a,b,i,j) = l2c(b,a,j,i)
                                 end do
                              end do
                           end do
                        end do

                        
                end subroutine update_L


                subroutine LH_1A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        integer, intent(in) :: iroot
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            l1a(sys%Nunocc_a,sys%Nocc_a), l1b(sys%Nunocc_b,sys%Nocc_b), &
                                            l2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            l2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            l2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: LH(sys%Nunocc_a,sys%Nocc_a)
                        real, allocatable :: Z(:,:), I1(:,:), I2(:,:), I3(:,:), I4(:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)
                        
                                if (iroot == 0) then
                                        LH = transpose(H1A%ou)
                                else
                                        LH = 0.0
                                end if

                                allocate(Z(nua,noa))
                                call einsum('ea,ei->ai',H1A%uu,l1a,Z)
                                LH = LH + Z
                                call einsum('im,am->ai',H1A%oo,l1a,Z)
                                LH = LH - Z
                                call einsum('eima,em->ai',H2A%uoou,l1a,Z)
                                LH = LH + Z
                                call einsum('ieam,em->ai',H2B%ouuo,l1b,Z)
                                LH = LH + Z
                                call einsum('fena,efin->ai',H2A%uuou,l2a,Z)
                                LH = LH + 0.5*Z
                                call einsum('efan,efin->ai',H2B%uuuo,l2b,Z)
                                LH = LH + Z
                                call einsum('finm,afmn->ai',H2A%uooo,l2a,Z)
                                LH = LH - 0.5*Z
                                call einsum('ifmn,afmn->ai',H2B%ouoo,l2b,Z)
                                LH = LH - Z

                                allocate(I1(nua,nua),I2(nua,nua),I3(noa,noa),I4(noa,noa))
                                call einsum('efmn,fgnm->ge',0.25*l2a,t2a,I1)
                                call einsum('efmn,egnm->gf',-0.25*l2a,t2a,I2)
                                call einsum('efmo,efno->mn',-0.25*l2a,t2a,I3)
                                call einsum('efmo,efnm->on',0.25*l2a,t2a,I4)
                                call einsum('ge,eiga->ai',I1,H2A%uouu,Z)
                                LH = LH + Z
                                call einsum('gf,figa->ai',I2,H2A%uouu,Z)
                                LH = LH + Z
                                call einsum('mn,nima->ai',I3,H2A%ooou,Z)
                                LH = LH + Z
                                call einsum('on,nioa->ai',I4,H2A%ooou,Z)
                                LH = LH + Z
                                deallocate(I1,I2,I3,I4)

                                allocate(I1(nob,nob),I2(nub,nub),I3(nua,nua),I4(noa,noa))
                                call einsum('abij,abin->jn',-1.0*l2b,t2b,I1)
                                call einsum('abij,afij->fb',l2b,t2b,I2)
                                call einsum('abij,fbij->fa',l2b,t2b,I3)
                                call einsum('abij,abnj->in',-1.0*l2b,t2b,I4)
                                call einsum('jn,mnej->em',I1,H2B%oouo,Z)
                                LH = LH + Z
                                call einsum('fb,mbef->em',I2,H2B%ouuu,Z)
                                LH = LH + Z
                                call einsum('fa,amfe->em',I3,H2A%uouu,Z)
                                LH = LH + Z
                                call einsum('in,nmie->em',I4,H2A%ooou,Z)
                                LH = LH + Z
                                deallocate(I1,I2,I3,I4)

                                allocate(I1(nub,nub),I2(nub,nub),I3(nob,nob),I4(nob,nob))
                                call einsum('abij,fbij->fa',0.25*l2c,t2c,I1)
                                call einsum('abij,faij->fb',-0.25*l2c,t2c,I2)
                                call einsum('abij,abnj->in',-0.25*l2c,t2c,I3)
                                call einsum('abij,abni->jn',0.25*l2c,t2c,I4)
                                call einsum('fa,maef->em',I1,H2B%ouuu,Z)
                                LH = LH + Z
                                call einsum('fb,mbef->em',I2,H2B%ouuu,Z)
                                LH = LH + Z
                                call einsum('in,mnei->em',I3,H2B%oouo,Z)
                                LH = LH + Z
                                call einsum('jn,mnej->em',I4,H2B%oouo,Z)
                                LH = LH + Z
                                deallocate(I1,I2,I3,I4)

                                deallocate(Z)
                        end associate


                end subroutine LH_1A

                subroutine LH_1B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        integer, intent(in) :: iroot
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            l1a(sys%Nunocc_a,sys%Nocc_a), l1b(sys%Nunocc_b,sys%Nocc_b), &
                                            l2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            l2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            l2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: LH(sys%Nunocc_b,sys%Nocc_b)
                        real, allocatable :: Z(:,:), I1(:,:), I2(:,:), I3(:,:), I4(:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                if (iroot == 0) then
                                        LH = transpose(H1B%ou)
                                else
                                        LH = 0.0
                                end if

                                allocate(Z(nub,nob))
                                call einsum('ea,ei->ai',H1B%uu,l1b,Z)
                                LH = LH + Z
                                call einsum('im,am->ai',H1B%oo,l1b,Z)
                                LH = LH - Z
                                call einsum('eima,em->ai',H2C%uoou,l1b,Z)
                                LH = LH + Z
                                call einsum('eima,em->ai',H2B%uoou,l1a,Z)
                                LH = LH + Z
                                call einsum('fena,efin->ai',H2C%uuou,l2c,Z)
                                LH = LH + 0.5*Z
                                call einsum('fena,feni->ai',H2B%uuou,l2b,Z)
                                LH = LH + Z
                                call einsum('finm,afmn->ai',H2C%uooo,l2c,Z)
                                LH = LH - 0.5*Z
                                call einsum('finm,fanm->ai',H2B%uooo,l2b,Z)
                                LH = LH - Z

                                allocate(I1(nub,nub),I2(nub,nub),I3(nob,nob),I4(nob,nob))
                                call einsum('efmn,fgnm->ge',0.25*l2c,t2c,I1)
                                call einsum('efmn,egnm->gf',-0.25*l2c,t2c,I2)
                                call einsum('efmo,efno->mn',-0.25*l2c,t2c,I3)
                                call einsum('efmo,efnm->on',0.25*l2c,t2c,I4)
                                call einsum('ge,eiga->ai',I1,H2C%uouu,Z)
                                LH = LH + Z
                                call einsum('gf,figa->ai',I2,H2C%uouu,Z)
                                LH = LH + Z
                                call einsum('mn,nima->ai',I3,H2C%ooou,Z)
                                LH = LH + Z
                                call einsum('on,nioa->ai',I4,H2C%ooou,Z)
                                LH = LH + Z
                                deallocate(I1,I2,I3,I4)

                                allocate(I1(nua,nua),I2(nub,nub),I3(noa,noa),I4(nob,nob))
                                call einsum('efmn,gfmn->ge',l2b,t2b,I1)
                                call einsum('fenm,fgnm->ge',l2b,t2b,I2)
                                call einsum('efmn,efon->mo',-1.0*l2b,t2b,I3)
                                call einsum('fenm,feno->mo',-1.0*l2b,t2b,I4)
                                call einsum('ge,eiga->ai',I1,H2B%uouu,Z)
                                LH = LH + Z
                                call einsum('ge,eiga->ai',I2,H2C%uouu,Z)
                                LH = LH + Z
                                call einsum('mo,oima->ai',I3,H2B%ooou,Z)
                                LH = LH + Z
                                call einsum('mo,oima->ai',I4,H2C%ooou,Z)
                                LH = LH + Z
                                deallocate(I1,I2,I3,I4)

                                allocate(I1(nua,nua),I2(nua,nua),I3(noa,noa),I4(noa,noa))
                                call einsum('abij,fbij->fa',0.25*l2a,t2a,I1)
                                call einsum('abij,faij->fb',-0.25*l2a,t2a,I2)
                                call einsum('abij,abnj->in',-0.25*l2a,t2a,I3)
                                call einsum('abij,abni->jn',0.25*l2a,t2a,I4)
                                call einsum('fa,amfe->em',I1,H2B%uouu,Z)
                                LH = LH + Z
                                call einsum('fb,bmfe->em',I2,H2B%uouu,Z)
                                LH = LH + Z
                                call einsum('in,nmie->em',I3,H2B%ooou,Z)
                                LH = LH + Z
                                call einsum('jn,nmje->em',I4,H2B%ooou,Z)
                                LH = LH + Z
                                deallocate(I1,I2,I3,I4)

                                deallocate(Z)
                        end associate

                end subroutine LH_1B

                subroutine LH_2A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        integer, intent(in) :: iroot
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            l1a(sys%Nunocc_a,sys%Nocc_a), l1b(sys%Nunocc_b,sys%Nocc_b), &
                                            l2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            l2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            l2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: LH(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a)
                        real, allocatable :: I1(:,:), I2(:,:), Q(:,:,:,:), Q1(:,:,:,:), Q2(:,:,:,:), Q3(:,:,:,:)
                        integer :: a, b, i, j

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                if (iroot == 0) then
                                        call reorder_stripe(4,shape(vA%oouu),size(vA%oouu),'3412',vA%oouu,LH)
                                else
                                        LH = 0.0
                                end if

                                do i = 1,noa
                                   do j = i+1,noa
                                      do a = 1,nua
                                         do b = a+1,nua
                                            LH(b,a,j,i) = LH(b,a,j,i) +H1A%ou(j,b)*l1a(a,i) &
                                                                      -H1A%ou(i,b)*l1a(a,j) &
                                                                      -H1A%ou(j,a)*l1a(b,i) &
                                                                      +H1A%ou(i,a)*l1a(b,j)
                                            LH(a,b,j,i) = -LH(b,a,j,i)
                                            LH(b,a,i,j) = -LH(b,a,j,i)
                                            LH(a,b,i,j) = LH(b,a,j,i)
                                         end do
                                      end do
                                   end do
                                end do


                                allocate(Q(nua,nua,noa,noa),Q1(nua,nua,noa,noa))

                                call einsum('ea,ebij->abij',H1A%uu,l2a,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                LH = LH + Q - Q1
                                call einsum('im,abmj->abij',H1A%oo,l2a,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                LH = LH - Q + Q1

                                allocate(I1(nua,nua),I2(nua,nua))
                                call einsum('afmn,efmn->ea',l2a,t2a,I1)
                                call einsum('bfmn,efmn->eb',l2a,t2a,I2)
                                call einsum('ea,ijeb->abij',I1,vA%oouu,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                LH = LH - 0.5*Q + 0.5*Q1
                                deallocate(I1,I2)

                                allocate(I1(nua,nua),I2(nua,nua))
                                call einsum('afmn,efmn->ea',l2b,t2b,I1)
                                call einsum('bfmn,efmn->eb',l2b,t2b,I2)
                                call einsum('ea,ijeb->abij',I1,vA%oouu,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                LH = LH - Q + Q1
                                deallocate(I1,I2)
                                
                                allocate(I1(noa,noa),I2(noa,noa))
                                call einsum('efin,efmn->im',l2a,t2a,I1)
                                call einsum('efjn,efmn->jm',l2a,t2a,I2)
                                call einsum('im,mjab->abij',I1,vA%oouu,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                LH = LH - 0.5*Q + 0.5*Q1
                                deallocate(I1,I2)

                                allocate(I1(noa,noa),I2(noa,noa))
                                call einsum('efin,efmn->im',l2b,t2b,I1)
                                call einsum('efjn,efmn->jm',l2b,t2b,I2)
                                call einsum('im,mjab->abij',I1,vA%oouu,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                LH = LH - Q + Q1
                                deallocate(I1,I2)

                                allocate(Q2(nua,nua,noa,noa),Q3(nua,nua,noa,noa))
                                call einsum('eima,ebmj->abij',H2A%uoou,l2a,Q)
                                call einsum('ieam,bejm->abij',H2B%ouuo,l2b,Q1)
                                Q = Q + Q1
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q2)
                                call reorder_stripe(4,shape(Q),size(Q),'2143',Q,Q3)
                                LH = LH + Q - Q1 - Q2 + Q3
                                deallocate(Q2,Q3)

                                call einsum('ijmn,abmn->abij',H2A%oooo,l2a,Q)
                                LH = LH + 0.5*Q
                                call einsum('efab,efij->abij',H2A%uuuu,l2a,Q)
                                LH = LH + 0.5*Q
                                call einsum('ejab,ei->abij',H2A%uouu,l1a,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                LH = LH + Q - Q1
                                call einsum('ijmb,am->abij',H2A%ooou,l1a,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                LH = LH - Q + Q1

                                deallocate(Q,Q1)


                        end associate


                end subroutine LH_2A
                
                subroutine LH_2B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        integer, intent(in) :: iroot
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            l1a(sys%Nunocc_a,sys%Nocc_a), l1b(sys%Nunocc_b,sys%Nocc_b), &
                                            l2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            l2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            l2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: LH(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b)
                        integer :: a, b, i, j
                        real, allocatable :: Q(:,:,:,:), I1(:,:), I2(:,:), I3(:,:), I4(:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                if (iroot == 0) then
                                        call reorder_stripe(4,shape(vB%oouu),size(vB%oouu),'3412',vB%oouu,LH)
                                else
                                        LH = 0.0
                                end if

                                do j = 1,nob
                                   do i = 1,noa
                                      do b = 1,nub
                                         do a = 1,nua
                                            LH(a,b,i,j) = LH(a,b,i,j) +H1B%ou(j,b)*l1a(a,i) &
                                                                      +H1A%ou(i,a)*l1b(b,j) 
                                         end do
                                      end do
                                   end do
                                end do


                                allocate(Q(nua,nub,noa,nob))

                                call einsum('ijmb,am->abij',H2B%ooou,l1a,Q)
                                LH = LH - Q
                                call einsum('ijam,bm->abij',H2B%oouo,l1b,Q)
                                LH = LH - Q
                                call einsum('ejab,ei->abij',H2B%uouu,l1a,Q)
                                LH = LH + Q
                                call einsum('ieab,ej->abij',H2B%ouuu,l1b,Q)
                                LH = LH + Q
                                call einsum('ijmn,abmn->abij',H2B%oooo,l2b,Q)
                                LH = LH + Q
                                call einsum('efab,efij->abij',H2B%uuuu,l2b,Q)
                                LH = LH + Q
                                call einsum('ejmb,aeim->abij',H2B%uoou,l2a,Q)
                                LH = LH + Q
                                call einsum('eima,ebmj->abij',H2A%uoou,l2b,Q)
                                LH = LH + Q
                                call einsum('ejmb,aeim->abij',H2C%uoou,l2b,Q)
                                LH = LH + Q
                                call einsum('ieam,ebmj->abij',H2B%ouuo,l2c,Q)
                                LH = LH + Q
                                call einsum('iemb,aemj->abij',H2B%ouou,l2b,Q)
                                LH = LH - Q
                                call einsum('ejam,ebim->abij',H2B%uouo,l2b,Q)
                                LH = LH - Q

                                allocate(I1(nua,nua),I2(nua,nua),I3(nub,nub),I4(nub,nub))
                                call einsum('abij,fbij->fa',l2a,t2a,I1)
                                call einsum('afmn,efmn->ea',l2b,t2b,I2)
                                call einsum('fbnm,fenm->eb',l2b,t2b,I3)
                                call einsum('bfmn,efmn->eb',l2c,t2c,I4)
                                call einsum('fa,nmfe->aenm',I1,vB%oouu,Q)
                                LH = LH - 0.5*Q
                                call einsum('ea,ijeb->abij',I2,vB%oouu,Q)
                                LH = LH - Q
                                call einsum('eb,ijae->abij',I3,vB%oouu,Q)
                                LH = LH - Q
                                call einsum('eb,ijae->abij',I4,vB%oouu,Q)
                                LH = LH - 0.5*Q
                                deallocate(I1,I2,I3,I4)

                                allocate(I1(noa,noa),I2(noa,noa),I3(nob,nob),I4(nob,nob))
                                call einsum('efin,efmn->im',l2a,t2a,I1)
                                call einsum('efin,efmn->im',l2b,t2b,I2)
                                call einsum('fenj,fenm->jm',l2b,t2b,I3)
                                call einsum('efjn,efmn->jm',l2c,t2c,I4)
                                call einsum('im,mjab->abij',I1,vB%oouu,Q)
                                LH = LH - 0.5*Q
                                call einsum('im,mjab->abij',I2,vB%oouu,Q)
                                LH = LH - Q
                                call einsum('jm,imab->abij',I3,vB%oouu,Q)
                                LH = LH - Q
                                call einsum('jm,imab->abij',I4,vB%oouu,Q)
                                LH = LH - 0.5*Q
                                deallocate(I1,I2,I3,I4)

                                call einsum('ea,ebij->abij',H1A%uu,l2b,Q)
                                LH = LH + Q
                                call einsum('eb,aeij->abij',H1B%uu,l2b,Q)
                                LH = LH + Q
                                call einsum('im,abmj->abij',H1A%oo,l2b,Q)
                                LH = LH - Q
                                call einsum('jm,abim->abij',H1B%oo,l2b,Q)
                                LH = LH - Q

                                deallocate(Q)

                        end associate

                end subroutine LH_2B

                subroutine LH_2C(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        integer, intent(in) :: iroot
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                            l1a(sys%Nunocc_a,sys%Nocc_a), l1b(sys%Nunocc_b,sys%Nocc_b), &
                                            l2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            l2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            l2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: LH(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, allocatable :: I1(:,:), I2(:,:), Q(:,:,:,:), Q1(:,:,:,:), Q2(:,:,:,:), Q3(:,:,:,:)
                        integer :: a, b, i, j

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                if (iroot == 0) then
                                        call reorder_stripe(4,shape(vC%oouu),size(vC%oouu),'3412',vC%oouu,LH)
                                else
                                        LH = 0.0
                                end if

                                do i = 1,nob
                                   do j = i+1,nob
                                      do a = 1,nub
                                         do b = a+1,nub
                                            LH(b,a,j,i) = LH(b,a,j,i) +H1B%ou(j,b)*l1b(a,i) &
                                                                      -H1B%ou(i,b)*l1b(a,j) &
                                                                      -H1B%ou(j,a)*l1b(b,i) &
                                                                      +H1B%ou(i,a)*l1b(b,j)
                                            LH(a,b,j,i) = -LH(b,a,j,i)
                                            LH(b,a,i,j) = -LH(b,a,j,i)
                                            LH(a,b,i,j) = LH(b,a,j,i)
                                         end do
                                      end do
                                   end do
                                end do


                                allocate(Q(nub,nub,nob,nob),Q1(nub,nub,nob,nob))

                                call einsum('ea,ebij->abij',H1B%uu,l2c,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                LH = LH + Q - Q1
                                call einsum('im,abmj->abij',H1B%oo,l2c,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                LH = LH - Q + Q1

                                allocate(I1(nub,nub),I2(nub,nub))
                                call einsum('afmn,efmn->ea',l2c,t2c,I1)
                                call einsum('bfmn,efmn->eb',l2c,t2c,I2)
                                call einsum('ea,ijeb->abij',I1,vC%oouu,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                LH = LH - 0.5*Q + 0.5*Q1
                                deallocate(I1,I2)

                                allocate(I1(nub,nub),I2(nub,nub))
                                call einsum('fanm,fenm->ea',l2b,t2b,I1)
                                call einsum('fbnm,fenm->eb',l2b,t2b,I2)
                                call einsum('ea,ijeb->abij',I1,vC%oouu,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                LH = LH - Q + Q1
                                deallocate(I1,I2)
                                
                                allocate(I1(nob,nob),I2(nob,nob))
                                call einsum('efin,efmn->im',l2c,t2c,I1)
                                call einsum('efjn,efmn->jm',l2c,t2c,I2)
                                call einsum('im,mjab->abij',I1,vC%oouu,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                LH = LH - 0.5*Q + 0.5*Q1
                                deallocate(I1,I2)

                                allocate(I1(nob,nob),I2(nob,nob))
                                call einsum('feni,fenm->im',l2b,t2b,I1)
                                call einsum('fenj,fenm->jm',l2b,t2b,I2)
                                call einsum('im,mjab->abij',I1,vC%oouu,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                LH = LH - Q + Q1
                                deallocate(I1,I2)

                                allocate(Q2(nub,nub,nob,nob),Q3(nub,nub,nob,nob))
                                call einsum('eima,ebmj->abij',H2C%uoou,l2c,Q)
                                call einsum('eima,ebmj->abij',H2B%uoou,l2b,Q1)
                                Q = Q + Q1
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q2)
                                call reorder_stripe(4,shape(Q),size(Q),'2143',Q,Q3)
                                LH = LH + Q - Q1 - Q2 + Q3
                                deallocate(Q2,Q3)

                                call einsum('ijmn,abmn->abij',H2C%oooo,l2c,Q)
                                LH = LH + 0.5*Q
                                call einsum('efab,efij->abij',H2C%uuuu,l2c,Q)
                                LH = LH + 0.5*Q
                                call einsum('ejab,ei->abij',H2C%uouu,l1b,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                LH = LH + Q - Q1
                                call einsum('ijmb,am->abij',H2C%ooou,l1b,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                LH = LH - Q + Q1

                                deallocate(Q,Q1)


                        end associate


                end subroutine LH_2C



end module leftccsd
