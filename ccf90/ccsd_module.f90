module ccsd_module

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use einsum_module, only: einsum
        use permutils, only: reorder_stripe

        implicit none

        contains

                subroutine ccsd(sys,fA,fB,vA,vB,vC,ndiis,maxit,shift,tol,t1a,t1b,t2a,t2b,t2c,Ecorr)

                        use cc_energy, only: calc_cc_energy
                        use diis, only: diis_xtrap

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        integer, intent(in) :: ndiis, maxit
                        real, intent(in) :: shift, tol
                        real, intent(out) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b, sys%Nocc_b), &
                                             t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                             t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                             t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                                             Ecorr
                        type(e1int_t) :: H1A, H1B
                        type(e2int_t) :: H2A, H2B, H2C
                        integer :: ndim, n1a, n1b, n2a, n2b, n2c, it_macro, it_micro
                        real :: Ecorr_old, dEcorr, resid
                        real, allocatable :: t(:), t_old(:), t_resid(:), t_list(:,:), t_resid_list(:,:)

                        n1a = sys%Nocc_a * sys%Nunocc_a
                        n1b = sys%Nocc_b * sys%Nunocc_b
                        n2a = sys%Nocc_a**2 * sys%Nunocc_a**2
                        n2b = sys%Nocc_a*sys%Nocc_b*sys%Nunocc_a*sys%Nunocc_b
                        n2c = sys%Nocc_b**2 * sys%Nunocc_b**2
                        ndim = n1a+n1b+n2a+n2b+n2c

                        allocate(t(ndim),t_old(ndim),t_resid(ndim),t_list(ndim,ndiis),t_resid_list(ndim,ndiis))

                        write(*,fmt=*) ''
                        write(*,fmt=*) 'CCSD ROUTINE'
                        write(*,fmt=*) ''

                        t = 0.0
                        Ecorr_old = 0.0
                        it_macro = 0
                        
                        write(*,fmt=*) '               Iteration      Ecorr                       Residuum'         
                        write(*,fmt=*) '------------------------------------------------------------------------'
                        do it_micro = 0,maxit

                                t_old = t
                                
                                t1a = reshape(t(1:n1a),(/sys%Nunocc_a,sys%Nocc_a/))
                                t1b = reshape(t(n1a+1:n1a+n1b),(/sys%Nunocc_b,sys%Nocc_b/))
                                t2a = reshape(t(n1a+n1b+1:n1a+n1b+n2a),(/sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a/))
                                t2b = reshape(t(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b),(/sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b/))
                                t2c = reshape(t(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c), &
                                        (/sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b/))  

                                call calc_cc_energy(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Ecorr)

                                call update_t1a(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,shift)
                                call update_t1b(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,shift)

                                call hbar_ccs_intermediates(sys,fA,fB,vA,vB,vC,t1a,t1b,H1A,H1B,H2A,H2B,H2C)

                                call update_t2a(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)
                                call update_t2b(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)
                                call update_t2c(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                                t(1:n1a) = reshape(t1a,(/n1a/))
                                t(n1a+1:n1a+n1b) = reshape(t1b,(/n1b/))
                                t(n1a+n1b+1:n1a+n1b+n2a) = reshape(t2a,(/n2a/))
                                t(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b) = reshape(t2b,(/n2b/))
                                t(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c) = reshape(t2c,(/n2c/))

                                t_resid = t - t_old
                                dEcorr = Ecorr - Ecorr_old
                                resid = norm2(t_resid)

                                if ( (resid <= tol) .AND. (abs(dEcorr) <= tol) ) then
                                    write(*,fmt=*) 'CCSD converged!'
                                    write(*,fmt=*) ''
                                    write(*,fmt=*) 'CCSD CORRELATION ENERGY (HARTREE) = ',Ecorr
                                    write(*,fmt=*) 'CCSD ENERGY (HARTREE) = ',Ecorr+sys%Escf
                                    exit
                                end if

                                t_list(:,mod(it_micro,ndiis)+1) = t
                                t_resid_list(:,mod(it_micro,ndiis)+1) = t_resid

                                if (mod(it_micro+1,ndiis) == 0) then
                                   it_macro = it_macro + 1
                                   write(*,fmt='(1x,a15,1x,i0)') 'DIIS cycle - ',it_macro
                                   call diis_xtrap(t,t_list,t_resid_list)
                                end if
                                    
                                write(*,fmt=*) it_micro,' ',Ecorr,' ',resid

                                Ecorr_old = Ecorr

                        end do

                        deallocate(t,t_old,t_resid,t_list,t_resid_list)

                end subroutine ccsd

                subroutine update_t1a(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,shift)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: t1b(sys%Nunocc_b, sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(in) :: shift
                        real, intent(inout) :: t1a(sys%Nunocc_a,sys%Nocc_a)
                        real :: denom, val
                        integer :: a, i
                        real :: I1A_vv(sys%Nunocc_a,sys%Nunocc_a), I1A_oo(sys%Nocc_a,sys%Nocc_a), &
                                h1A_ov(sys%Nocc_a,sys%Nunocc_a), h1B_ov(sys%Nocc_b,sys%Nunocc_b), &
                                h1A_oo(sys%Nocc_a,sys%Nocc_a), &
                                h2A_ooov(sys%Nocc_a,sys%Nocc_a,sys%Nocc_a,sys%Nunocc_a), &
                                h2B_ooov(sys%Nocc_a,sys%Nocc_b,sys%Nocc_a,sys%Nunocc_b), &
                                h2A_vovv(sys%Nunocc_a,sys%Nocc_a,sys%Nunocc_a,sys%Nunocc_a), &
                                h2B_vovv(sys%Nunocc_a,sys%Nocc_b,sys%Nunocc_a,sys%Nunocc_b), &
                                X1A(sys%Nunocc_a,sys%Nocc_a)
                        real, allocatable :: Z1(:,:), Q1(:,:,:,:)

                        associate(noa=>sys%Nocc_a,nob=>sys%Nocc_b,nua=>sys%Nunocc_a,nub=>sys%Nunocc_b)

                                I1A_vv = fA%uu
                                allocate(Z1(nua,nua))
                                call einsum('anef,fn->ae',vA%uouu,t1a,Z1)
                                I1A_vv = I1A_vv + Z1
                                call einsum('anef,fn->ae',vB%uouu,t1b,Z1)
                                I1A_vv = I1A_vv + Z1
                                deallocate(Z1)

                                I1A_oo = fA%oo
                                allocate(Z1(noa,noa))
                                call einsum('mnif,fn->mi',vA%ooou,t1a,Z1)
                                I1A_oo = I1A_oo + Z1
                                call einsum('mnif,fn->mi',vB%ooou,t1b,Z1)
                                I1A_oo = I1A_oo + Z1
                                deallocate(Z1)

                                h1A_ov = fA%ou
                                allocate(Z1(noa,nua))
                                call einsum('mnef,fn->me',vA%oouu,t1a,Z1)
                                h1A_ov = h1A_ov + Z1
                                call einsum('mnef,fn->me',vB%oouu,t1b,Z1)
                                h1A_ov = h1A_ov + Z1
                                deallocate(Z1)

                                h1B_ov = fB%ou
                                allocate(Z1(nob,nub))
                                call einsum('nmfe,fn->me',vB%oouu,t1a,Z1)
                                h1B_ov = h1B_ov + Z1
                                call einsum('mnef,fn->me',vC%oouu,t1b,Z1)
                                h1B_ov = h1B_ov + Z1
                                deallocate(Z1)

                                h1A_oo = I1A_oo
                                allocate(Z1(noa,noa))
                                call einsum('me,ei->mi',h1A_ov,t1a,Z1)
                                h1A_oo = h1A_oo + Z1
                                deallocate(Z1)

                                h2A_ooov = vA%ooou
                                allocate(Q1(noa,noa,noa,nua))
                                call einsum('mnfe,fi->mnie',vA%oouu,t1a,Q1)
                                h2A_ooov = h2A_ooov + Q1
                                deallocate(Q1)

                                h2B_ooov = vB%ooou
                                allocate(Q1(noa,nob,noa,nub))
                                call einsum('mnfe,fi->mnie',vB%oouu,t1a,Q1)
                                h2B_ooov = h2B_ooov + Q1
                                deallocate(Q1)

                                h2A_vovv = vA%uouu
                                allocate(Q1(nua,noa,nua,nua))
                                call einsum('mnfe,an->amef',vA%oouu,t1a,Q1)
                                h2A_vovv = h2A_vovv - Q1
                                deallocate(Q1)

                                h2B_vovv = vB%uouu
                                allocate(Q1(nua,nob,nua,nub))
                                call einsum('nmef,an->amef',vB%oouu,t1a,Q1)
                                h2B_vovv = h2B_vovv - Q1
                                deallocate(Q1)

                                X1A = fA%uo
                                allocate(Z1(nua,noa))
                                call einsum('mi,am->ai',h1A_oo,t1a,Z1)
                                X1A = X1A - Z1
                                call einsum('ae,ei->ai',I1A_vv,t1a,Z1)
                                X1A = X1A + Z1
                                call einsum('anif,fn->ai',vA%uoou,t1a,Z1)
                                X1A = X1A + Z1
                                call einsum('anif,fn->ai',vB%uoou,t1b,Z1)
                                X1A = X1A + Z1
                                call einsum('me,aeim->ai',h1A_ov,t2a,Z1)
                                X1A = X1A + Z1
                                call einsum('me,aeim->ai',h1B_ov,t2b,Z1)
                                X1A = X1A + Z1
                                call einsum('mnif,afmn->ai',h2A_ooov,t2a,Z1)
                                X1A = X1A - 0.5*Z1
                                call einsum('mnif,afmn->ai',h2B_ooov,t2b,Z1)
                                X1A = X1A - Z1
                                call einsum('anef,efin->ai',h2A_vovv,t2a,Z1)
                                X1A = X1A + 0.5*Z1
                                call einsum('anef,efin->ai',h2B_vovv,t2b,Z1)
                                X1A = X1A + Z1
                                deallocate(Z1)

                        end associate

                        do i = 1,sys%Nocc_a
                           do a = 1,sys%Nunocc_a
                              denom = fA%oo(i,i) - fA%uu(a,a)
                              val = X1A(a,i)
                              t1a(a,i) = t1a(a,i) + val/(denom-shift)
                           end do
                        end do

                end subroutine update_t1a


                subroutine update_t1b(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,shift)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: t1a(sys%Nunocc_a, sys%Nocc_a), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(in) :: shift
                        real, intent(inout) :: t1b(sys%Nunocc_b,sys%Nocc_a)
                        real :: denom, val
                        integer :: a, i
                        real :: I1B_vv(sys%Nunocc_b,sys%Nunocc_b), I1B_oo(sys%Nocc_b,sys%Nocc_b), &
                                h1A_ov(sys%Nocc_a,sys%Nunocc_a), h1B_ov(sys%Nocc_b,sys%Nunocc_b), &
                                h1B_oo(sys%Nocc_b,sys%Nocc_b), &
                                h2C_ooov(sys%Nocc_b,sys%Nocc_b,sys%Nocc_b,sys%Nunocc_b), &
                                h2B_oovo(sys%Nocc_a,sys%Nocc_b,sys%Nunocc_a,sys%Nocc_b), &
                                h2C_vovv(sys%Nunocc_b,sys%Nocc_b,sys%Nunocc_b,sys%Nunocc_b), &
                                h2B_ovvv(sys%Nocc_a,sys%Nunocc_b,sys%Nunocc_a,sys%Nunocc_b), &
                                X1B(sys%Nunocc_b,sys%Nocc_b)
                        real, allocatable :: Z1(:,:), Q1(:,:,:,:)

                        associate(noa=>sys%Nocc_a,nob=>sys%Nocc_b,nua=>sys%Nunocc_a,nub=>sys%Nunocc_b)
                                
                                I1B_vv = fB%uu
                                allocate(Z1(nub,nub))
                                call einsum('anef,fn->ae',vC%uouu,t1b,Z1)
                                I1B_vv = I1B_vv + Z1
                                call einsum('nafe,fn->ae',vB%ouuu,t1a,Z1)
                                I1B_vv = I1B_vv + Z1
                                deallocate(Z1)

                                I1B_oo = fB%oo
                                allocate(Z1(nob,nob))
                                call einsum('mnif,fn->mi',vC%ooou,t1b,Z1)
                                I1B_oo = I1B_oo + Z1
                                call einsum('nmfi,fn->mi',vB%oouo,t1a,Z1)
                                I1B_oo = I1B_oo + Z1
                                deallocate(Z1)

                                h1A_ov = fA%ou
                                allocate(Z1(noa,nua))
                                call einsum('mnef,fn->me',vA%oouu,t1a,Z1)
                                h1A_ov = h1A_ov + Z1
                                call einsum('mnef,fn->me',vB%oouu,t1b,Z1)
                                h1A_ov = h1A_ov + Z1
                                deallocate(Z1)

                                h1B_ov = fB%ou
                                allocate(Z1(nob,nub))
                                call einsum('nmfe,fn->me',vB%oouu,t1a,Z1)
                                h1B_ov = h1B_ov + Z1
                                call einsum('mnef,fn->me',vC%oouu,t1b,Z1)
                                h1B_ov = h1B_ov + Z1
                                deallocate(Z1)

                                h1B_oo = I1B_oo
                                allocate(Z1(nob,nob))
                                call einsum('me,ei->mi',h1B_ov,t1b,Z1)
                                h1B_oo = h1B_oo + Z1
                                deallocate(Z1)

                                h2C_ooov = vC%ooou
                                allocate(Q1(nob,nob,nob,nub))
                                call einsum('mnfe,fi->mnie',vC%oouu,t1b,Q1)
                                h2C_ooov = h2C_ooov + Q1
                                deallocate(Q1)

                                h2B_oovo = vB%oouo
                                allocate(Q1(noa,nob,nua,nob))
                                call einsum('nmef,fi->nmei',vB%oouu,t1b,Q1)
                                h2B_oovo = h2B_oovo + Q1
                                deallocate(Q1)

                                h2C_vovv = vC%uouu
                                allocate(Q1(nub,nob,nub,nub))
                                call einsum('mnfe,an->amef',vC%oouu,t1b,Q1)
                                h2C_vovv = h2C_vovv - Q1
                                deallocate(Q1)

                                h2B_ovvv = vB%ouuu
                                allocate(Q1(noa,nub,nua,nub))
                                call einsum('mnfe,an->mafe',vB%oouu,t1b,Q1)
                                h2B_ovvv = h2B_ovvv - Q1
                                deallocate(Q1)

                                X1B = fB%uo
                                allocate(Z1(nub,nob))
                                call einsum('mi,am->ai',h1B_oo,t1b,Z1)
                                X1B = X1B - Z1
                                call einsum('ae,ei->ai',I1B_vv,t1b,Z1)
                                X1B = X1B + Z1
                                call einsum('anif,fn->ai',vC%uoou,t1b,Z1)
                                X1B = X1B + Z1
                                call einsum('nafi,fn->ai',vB%ouuo,t1a,Z1)
                                X1B = X1B + Z1
                                call einsum('me,eami->ai',h1A_ov,t2b,Z1)
                                X1B = X1B + Z1
                                call einsum('me,aeim->ai',h1B_ov,t2c,Z1)
                                X1B = X1B + Z1
                                call einsum('mnif,afmn->ai',h2C_ooov,t2c,Z1)
                                X1B = X1B - 0.5*Z1
                                call einsum('nmfi,fanm->ai',h2B_oovo,t2b,Z1)
                                X1B = X1B - Z1
                                call einsum('anef,efin->ai',h2C_vovv,t2c,Z1)
                                X1B = X1B + 0.5*Z1
                                call einsum('nafe,feni->ai',h2B_ovvv,t2b,Z1)
                                X1B = X1B + Z1
                                deallocate(Z1)

                        end associate

                        do i = 1,sys%Nocc_b
                           do a = 1,sys%Nunocc_b
                              denom = fB%oo(i,i) - fB%uu(a,a)
                              val = X1B(a,i)
                              t1b(a,i) = t1b(a,i) + val/(denom-shift)
                           end do
                        end do

                end subroutine update_t1b


                subroutine update_t2a(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(in) :: shift
                        integer :: a, b, i, j
                        real, intent(inout) :: t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a)
                        real :: I1A_oo(sys%Nocc_a,sys%Nocc_a), I1A_vv(sys%Nunocc_a,sys%Nunocc_a), &
                                I2A_voov(sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a,sys%Nunocc_a), &
                                I2A_oooo(sys%Nocc_a,sys%Nocc_a,sys%Nocc_a,sys%Nocc_a), &
                                I2B_voov(sys%Nunocc_a,sys%Nocc_b,sys%Nocc_a,sys%Nunocc_b), &
                                X2A(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                denom, val
                        real, allocatable :: Z1(:,:), Q1(:,:,:,:), &
                                             D1(:,:,:,:), D2(:,:,:,:), D3(:,:,:,:), D4(:,:,:,:), &
                                             D5(:,:,:,:), D6(:,:,:,:)

                        associate(noa=>sys%Nocc_a,nob=>sys%Nocc_b,nua=>sys%Nunocc_a,nub=>sys%Nunocc_b)


                                I1A_oo = H1A%oo
                                allocate(Z1(noa,noa))
                                call einsum('mnef,efin->mi',vA%oouu,t2a,Z1)
                                I1A_oo = I1A_oo + 0.5*Z1
                                call einsum('mnef,efin->mi',vB%oouu,t2b,Z1)
                                I1A_oo = I1A_oo + Z1
                                deallocate(Z1)

                                I1A_vv = H1A%uu
                                allocate(Z1(nua,nua))
                                call einsum('mnef,afmn->ae',vA%oouu,t2a,Z1)
                                I1A_vv = I1A_vv - 0.5*Z1
                                call einsum('mnef,afmn->ae',vB%oouu,t2b,Z1)
                                I1A_vv = I1A_vv - Z1
                                deallocate(Z1)

                                I2A_voov = H2A%uoou
                                allocate(Q1(nua,noa,noa,nua))
                                call einsum('mnef,afin->amie',vA%oouu,t2a,Q1)
                                I2A_voov = I2A_voov + 0.5*Q1
                                call einsum('mnef,afin->amie',vB%oouu,t2b,Q1)
                                I2A_voov = I2A_voov + Q1
                                deallocate(Q1)

                                I2A_oooo = H2A%oooo
                                allocate(Q1(noa,noa,noa,noa))
                                call einsum('mnef,efij->mnij',vA%oouu,t2a,Q1)
                                I2A_oooo = I2A_oooo + 0.5*Q1
                                deallocate(Q1)

                                I2B_voov = H2B%uoou
                                !I2B_voov = 0.0
                                allocate(Q1(nua,nob,noa,nub))
                                call einsum('mnef,afin->amie',vC%oouu,t2b,Q1)
                                I2B_voov = I2B_voov + 0.5*Q1
                                deallocate(Q1)

                                X2A = vA%uuoo
                                allocate(D1(nua,nua,noa,noa),D2(nua,nua,noa,noa),D3(nua,nua,noa,noa), &
                                         D4(nua,nua,noa,noa),D5(nua,nua,noa,noa),D6(nua,nua,noa,noa))

                                call einsum('amij,bm->abij',-1.0*H2A%uooo,t1a,D1)
                                call einsum('abie,ej->abij',H2A%uuou,t1a,D2)
                                call einsum('ae,ebij->abij',I1A_vv,t2a,D3)
                                call einsum('mi,abmj->abij',-1.0*I1A_oo,t2a,D4)
                                call einsum('amie,ebmj->abij',I2A_voov,t2a,D5)
                                call einsum('amie,bejm->abij',I2B_voov,t2b,D6)

                                allocate(Q1(nua,nua,noa,noa))
                                call einsum('abef,efij->abij',H2A%uuuu,t2a,Q1)
                                X2A = X2A + 0.5*Q1
                                call einsum('mnij,abmn->abij',I2A_oooo,t2a,Q1)
                                X2A = X2A + 0.5*Q1
                                deallocate(Q1)

                                ! these diagrams are A(ab)
                                D1 = D1 + D3
                                ! these diagrams are A(ij)
                                D2 = D2 + D4
                                ! these diagrams are A(ab)A(ij)
                                D5 = D5 + D6

                                deallocate(D3,D4,D6)

                        end associate
                
                        do i = 1,sys%Nocc_a
                           do j = i+1,sys%Nocc_a
                              do a = 1,sys%Nunocc_a
                                 do b = a+1,sys%Nunocc_a
                                    
                                    denom = fA%oo(i,i)+fA%oo(j,j)-fA%uu(a,a)-fA%uu(b,b)
                                    val = X2A(b,a,j,i) &
                                          + D1(b,a,j,i) - D1(a,b,j,i) & ! A(ab)
                                          + D2(b,a,j,i) - D2(b,a,i,j) & ! A(ij)
                                          + D5(b,a,j,i) - D5(a,b,j,i) - D5(b,a,i,j) + D5(a,b,i,j) ! A(ab)A(ij)
                                    t2a(b,a,j,i) = t2a(b,a,j,i) + val/(denom-shift)
                                    t2a(a,b,j,i) = -t2a(b,a,j,i)
                                    t2a(b,a,i,j) = -t2a(b,a,j,i)
                                    t2a(a,b,i,j) = t2a(b,a,j,i)

                                 end do
                              end do
                           end do
                        end do

                        deallocate(D1,D2,D5)

                end subroutine update_t2a


                subroutine update_t2b(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(in) :: shift
                        real :: denom, val
                        integer :: a, b, i, j
                        real, intent(inout) :: t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b)
                        real :: I1A_vv(sys%Nunocc_a,sys%Nunocc_a), I1B_vv(sys%Nunocc_b,sys%Nunocc_b), &
                                I1A_oo(sys%Nocc_a,sys%Nocc_a), I1B_oo(sys%Nocc_b,sys%Nocc_b), &
                                I2A_voov(sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a,sys%Nunocc_a), &
                                I2B_voov(sys%Nunocc_a,sys%Nocc_b,sys%Nocc_a,sys%Nunocc_b), &
                                I2B_oooo(sys%Nocc_a,sys%Nocc_b,sys%Nocc_a,sys%Nocc_b), &
                                I2B_vovo(sys%Nunocc_a,sys%Nocc_b,sys%Nunocc_a,sys%Nocc_b), &
                                X2B(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b)
                        real, allocatable :: Z1(:,:), Q1(:,:,:,:)
               
                        associate(noa=>sys%Nocc_a,nob=>sys%Nocc_b,nua=>sys%Nunocc_a,nub=>sys%Nunocc_b)

                                I1A_vv = H1A%uu
                                allocate(Z1(nua,nua))
                                call einsum('mnef,afmn->ae',vA%oouu,t2a,Z1)
                                I1A_vv = I1A_vv - 0.5*Z1
                                call einsum('mnef,afmn->ae',vB%oouu,t2b,Z1)
                                I1A_vv = I1A_vv - Z1
                                deallocate(Z1)

                                I1B_vv = H1B%uu
                                allocate(Z1(nub,nub))
                                call einsum('nmfe,fbnm->be',vB%oouu,t2b,Z1)
                                I1B_vv = I1B_vv - Z1
                                call einsum('mnef,fbnm->be',vC%oouu,t2c,Z1)
                                I1B_vv = I1B_vv - 0.5*Z1
                                deallocate(Z1)

                                I1A_oo = H1A%oo
                                allocate(Z1(noa,noa))
                                call einsum('mnef,efin->mi',vA%oouu,t2a,Z1)
                                I1A_oo = I1A_oo + 0.5*Z1
                                call einsum('mnef,efin->mi',vB%oouu,t2b,Z1)
                                I1A_oo = I1A_oo + Z1
                                deallocate(Z1)

                                I1B_oo = H1B%oo
                                allocate(Z1(nob,nob))
                                call einsum('nmfe,fenj->mj',vB%oouu,t2b,Z1)
                                I1B_oo = I1B_oo + Z1
                                call einsum('mnef,efjn->mj',vC%oouu,t2c,Z1)
                                I1B_oo = I1B_oo + 0.5*Z1
                                deallocate(Z1)

                                I2A_voov = H2A%uoou
                                allocate(Q1(nua,noa,noa,nua))
                                call einsum('mnef,aeim->anif',vA%oouu,t2a,Q1)
                                I2A_voov = I2A_voov + Q1
                                call einsum('nmfe,aeim->anif',vB%oouu,t2b,Q1)
                                I2A_voov = I2A_voov + Q1
                                deallocate(Q1)

                                I2B_voov = H2B%uoou
                                allocate(Q1(nua,nob,noa,nub))
                                call einsum('mnef,aeim->anif',vB%oouu,t2a,Q1)
                                I2B_voov = I2B_voov + Q1
                                call einsum('mnef,aeim->anif',vC%oouu,t2b,Q1)
                                I2B_voov = I2B_voov + Q1
                                deallocate(Q1)
   
                                I2B_oooo = H2B%oooo
                                allocate(Q1(noa,nob,noa,nob))
                                call einsum('mnef,efij->mnij',vB%oouu,t2b,Q1)
                                I2B_oooo = I2B_oooo + Q1
                                deallocate(Q1)

                                I2B_vovo = H2B%uouo
                                allocate(Q1(nua,nob,nua,nob))
                                call einsum('mnef,afmj->anej',vB%oouu,t2b,Q1)
                                I2B_vovo = I2B_vovo - Q1
                                deallocate(Q1) 

                                X2B = vB%uuoo
                                allocate(Q1(nua,nub,noa,nob))
                                call einsum('mbij,am->abij',H2B%ouoo,t1a,Q1)
                                X2B = X2B - Q1
                                call einsum('amij,bm->abij',H2B%uooo,t1b,Q1)
                                X2B = X2B - Q1
                                call einsum('abej,ei->abij',H2B%uuuo,t1a,Q1)
                                X2B = X2B + Q1
                                call einsum('abie,ej->abij',H2B%uuou,t1b,Q1)
                                X2B = X2B + Q1
                                call einsum('ae,ebij->abij',I1A_vv,t2b,Q1)
                                X2B = X2B + Q1
                                call einsum('be,aeij->abij',I1B_vv,t2b,Q1)
                                X2B = X2B + Q1 
                                call einsum('mi,abmj->abij',I1A_oo,t2b,Q1)
                                X2B = X2B - Q1                
                                call einsum('mj,abim->abij',I1B_oo,t2b,Q1)
                                X2B = X2B - Q1    
                                call einsum('amie,ebmj->abij',I2A_voov,t2b,Q1)
                                X2B = X2B + Q1
                                call einsum('amie,ebmj->abij',I2B_voov,t2c,Q1)
                                X2B = X2B + Q1
                                call einsum('mbej,aeim->abij',H2B%ouuo,t2a,Q1)
                                X2B = X2B + Q1
                                call einsum('bmje,aeim->abij',H2C%uoou,t2b,Q1)
                                X2B = X2B + Q1        
                                call einsum('mbie,aemj->abij',H2B%ouou,t2b,Q1)
                                X2B = X2B - Q1
                                call einsum('amej,ebim->abij',I2B_vovo,t2b,Q1)
                                X2B = X2B - Q1     
                                call einsum('mnij,abmn->abij',I2B_oooo,t2b,Q1)
                                X2B = X2B + Q1     
                                call einsum('abef,efij->abij',H2B%uuuu,t2b,Q1)         
                                X2B = X2B + Q1
                                deallocate(Q1)

                        end associate
        
                        do j = 1,sys%Nocc_b
                           do i = 1,sys%Nocc_a
                              do b = 1,sys%Nunocc_b
                                 do a = 1,sys%Nunocc_a
                                    denom = fA%oo(i,i)+fB%oo(j,j)-fA%uu(a,a)-fB%uu(b,b)
                                    val = X2B(a,b,i,j)
                                    t2b(a,b,i,j) = t2b(a,b,i,j) + val/(denom-shift)
                                 end do
                              end do
                           end do
                        end do

                end subroutine update_t2b


                subroutine update_t2c(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b)
                        real, intent(in) :: shift
                        real, intent(inout) :: t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real :: denom, val
                        integer :: a, b, i, j
                        real :: I1B_oo(sys%Nocc_b,sys%Nocc_b), I1B_vv(sys%Nunocc_b,sys%Nunocc_b), &
                                I2C_voov(sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b,sys%Nunocc_b), &
                                I2C_oooo(sys%Nocc_b,sys%Nocc_b,sys%Nocc_b,sys%Nocc_b), &
                                I2B_ovvo(sys%Nocc_a,sys%Nunocc_b,sys%Nunocc_a,sys%Nocc_b), &
                                X2C(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, allocatable :: Z1(:,:), Q1(:,:,:,:), &
                                             D1(:,:,:,:), D2(:,:,:,:), D3(:,:,:,:), D4(:,:,:,:), &
                                             D5(:,:,:,:), D6(:,:,:,:)

                        associate(noa=>sys%Nocc_a,nob=>sys%Nocc_b,nua=>sys%Nunocc_a,nub=>sys%Nunocc_b)

                                I1B_oo = H1B%oo
                                allocate(Z1(nob,nob))
                                call einsum('mnef,efin->mi',vC%oouu,t2c,Z1)
                                I1B_oo = I1B_oo + 0.5*Z1
                                call einsum('nmfe,feni->mi',vB%oouu,t2b,Z1)
                                I1B_oo = I1B_oo + Z1
                                deallocate(Z1)

                                I1B_vv = H1B%uu
                                allocate(Z1(nub,nub))
                                call einsum('mnef,afmn->ae',vC%oouu,t2c,Z1)
                                I1B_vv = I1B_vv - 0.5*Z1
                                call einsum('nmfe,fanm->ae',vB%oouu,t2b,Z1)
                                I1B_vv = I1B_vv - Z1
                                deallocate(Z1)

                                I2C_oooo = H2C%oooo
                                allocate(Q1(nob,nob,nob,nob))
                                call einsum('mnef,efij->mnij',vC%oouu,t2c,Q1)
                                I2C_oooo = I2C_oooo + 0.5*Q1
                                deallocate(Q1)

                                I2B_ovvo = H2B%ouuo
                                allocate(Q1(noa,nub,nua,nob))
                                call einsum('mnef,afin->maei',vB%oouu,t2c,Q1)
                                I2B_ovvo = I2B_ovvo + Q1
                                call einsum('mnef,fani->maei',vA%oouu,t2b,Q1)
                                I2B_ovvo = I2B_ovvo + 0.5*Q1
                                deallocate(Q1)

                                I2C_voov = H2C%uoou
                                allocate(Q1(nub,nob,nob,nub))
                                call einsum('mnef,afin->amie',vC%oouu,t2c,Q1)
                                I2C_voov = I2C_voov + 0.5*Q1
                                deallocate(Q1) 

                                X2C = vC%uuoo
                                allocate(D1(nub,nub,nob,nob),D2(nub,nub,nob,nob),D3(nub,nub,nob,nob), &
                                         D4(nub,nub,nob,nob),D5(nub,nub,nob,nob),D6(nub,nub,nob,nob))

                                call einsum('bmji,am->abij',-1.0*H2C%uooo,t1b,D1)
                                call einsum('abej,ei->abij',H2C%uuuo,t1b,D2)
                                call einsum('ae,ebij->abij',I1B_vv,t2c,D3)
                                call einsum('mi,abmj->abij',-1.0*I1B_oo,t2c,D4)
                                call einsum('amie,ebmj->abij',I2C_voov,t2c,D5)
                                call einsum('maei,ebmj->abij',I2B_ovvo,t2b,D6)

                                allocate(Q1(nub,nub,nob,nob))
                                call einsum('abef,efij->abij',H2C%uuuu,t2c,Q1)
                                X2C = X2C + 0.5*Q1
                                call einsum('mnij,abmn->abij',I2C_oooo,t2c,Q1)
                                X2C = X2C + 0.5*Q1
                                deallocate(Q1)

                                D1 = D1 + D3 ! diagrams with A(ab)
                                D2 = D2 + D4 ! diagrams with A(ij)
                                D5 = D5 + D6 ! diagrams with A(ij)A(ab(

                                deallocate(D3,D4,D6)

                        end associate

                        do i = 1,sys%Nocc_b
                           do j = i+1,sys%Nocc_b
                              do a = 1,sys%Nunocc_b
                                 do b = a+1,sys%Nunocc_b
                                    
                                    denom = fB%oo(i,i)+fB%oo(j,j)-fB%uu(a,a)-fB%uu(b,b)
                                    val = X2C(b,a,j,i) &
                                          + D1(b,a,j,i) - D1(a,b,j,i) & ! A(ab)
                                          + D2(b,a,j,i) - D2(b,a,i,j) & ! A(ij)
                                          + D5(b,a,j,i) - D5(a,b,j,i) - D5(b,a,i,j) + D5(a,b,i,j) ! A(ab)A(ij)
                                    t2c(b,a,j,i) = t2c(b,a,j,i) + val/(denom-shift)
                                    t2c(a,b,j,i) = -t2c(b,a,j,i)
                                    t2c(b,a,i,j) = -t2c(b,a,j,i)
                                    t2c(a,b,i,j) = t2c(b,a,j,i)

                                 end do
                              end do
                           end do
                        end do

                        deallocate(D1,D2,D5)

                end subroutine update_t2c


                subroutine hbar_ccs_intermediates(sys,fA,fB,vA,vB,vC,t1a,t1b,H1A,H1B,H2A,H2B,H2C)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b, sys%Nocc_b)
                        type(e1int_t), intent(out) :: H1A, H1B
                        type(e2int_t), intent(out) :: H2A, H2B, H2C
                        real, allocatable :: Z(:,:), Q(:,:,:,:), Q2(:,:,:,:), Q3(:,:,:,:)

                        associate(noa=>sys%Nocc_a, nob=>sys%Nocc_b, nua=>sys%Nunocc_a, nub=>sys%Nunocc_b)

                                ! H1A(ou)
                                H1A%ou = fA%ou
                                allocate(Z(noa,nua))
                                call einsum('mnef,fn->me',vA%oouu,t1a,Z)
                                H1A%ou = H1A%ou + Z
                                call einsum('mnef,fn->me',vB%oouu,t1b,Z)
                                H1A%ou = H1A%ou + Z
                                deallocate(Z)

                                ! H1B(ou)
                                H1B%ou = fB%ou
                                allocate(Z(nob,nub))
                                call einsum('nmfe,fn->me',vB%oouu,t1a,Z)
                                H1B%ou = H1B%ou + Z
                                call einsum('mnef,fn->me',vC%oouu,t1b,Z)
                                H1B%ou = H1B%ou + Z
                                deallocate(Z)

                                ! H1A(oo)
                                H1A%oo = fA%oo
                                allocate(Z(noa,noa))
                                call einsum('mnif,fn->mi',vA%ooou,t1a,Z)
                                H1A%oo = H1A%oo + Z
                                call einsum('mnif,fn->mi',vB%ooou,t1b,Z)
                                H1A%oo = H1A%oo + Z
                                call einsum('me,ei->mi',H1A%ou,t1a,Z)
                                H1A%oo = H1A%oo + Z
                                deallocate(Z)

                                ! H1B(oo)
                                H1B%oo = fB%oo
                                allocate(Z(nob,nob))
                                call einsum('mnif,fn->mi',vC%ooou,t1b,Z)
                                H1B%oo = H1B%oo + Z
                                call einsum('nmfi,fn->mi',vB%oouo,t1a,Z)
                                H1B%oo = H1B%oo + Z
                                call einsum('me,ei->mi',H1B%ou,t1b,Z)
                                H1B%oo = H1B%oo + Z
                                deallocate(Z)

                                ! H1A(uu)
                                H1A%uu = fA%uu
                                allocate(Z(nua,nua))
                                call einsum('anef,fn->ae',vA%uouu,t1a,Z)
                                H1A%uu = H1A%uu + Z
                                call einsum('anef,fn->ae',vB%uouu,t1b,Z)
                                H1A%uu = H1A%uu + Z
                                call einsum('me,am->ae',H1A%ou,t1a,Z)
                                H1A%uu = H1A%uu - Z
                                deallocate(Z)

                                ! H1B(uu)
                                H1B%uu = fB%uu
                                allocate(Z(nub,nub))
                                call einsum('anef,fn->ae',vC%uouu,t1b,Z)
                                H1B%uu = H1B%uu + Z
                                call einsum('nafe,fn->ae',vB%ouuu,t1a,Z)
                                H1B%uu = H1B%uu + Z
                                call einsum('me,am->ae',H1B%ou,t1b,Z)
                                H1B%uu = H1B%uu - Z
                                deallocate(Z)

                                ! H2A(oooo)
                                H2A%oooo = vA%oooo
                                allocate(Q(noa,noa,noa,noa),Q2(noa,noa,nua,noa),Q3(noa,noa,noa,noa))
                                call einsum('mnef,fj->mnej',vA%oouu,t1a,Q2)
                                call einsum('mnej,ei->mnij',vA%oouo+0.5*Q2,t1a,Q)
                                H2A%oooo = H2A%oooo + Q
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q3)
                                H2A%oooo = H2A%oooo - Q3
                                deallocate(Q,Q2,Q3)

                                ! H2A(uuuu)
                                H2A%uuuu = vA%uuuu
                                allocate(Q(nua,nua,nua,nua),Q2(noa,nua,nua,nua),Q3(nua,nua,nua,nua))
                                call einsum('mnef,bn->mbef',vA%oouu,t1a,Q2)
                                call einsum('mbef,am->abef',vA%ouuu-0.5*Q2,t1a,Q)
                                H2A%uuuu = H2A%uuuu - Q
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q3)
                                H2A%uuuu = H2A%uuuu + Q3
                                deallocate(Q,Q2,Q3)

                                ! H2A(uooo)
                                H2A%uooo = vA%uooo
                                allocate(Q(nua,noa,noa,noa),Q2(nua,noa,noa,nua),Q3(nua,noa,nua,noa))
                                call einsum('nmij,an->amij',vA%oooo,t1a,Q)
                                H2A%uooo = H2A%uooo - 0.5*Q
                                call einsum('amef,ei->amif',vA%uouu,t1a,Q2)
                                call einsum('amif,fj->amij',Q2,t1a,Q)
                                H2A%uooo = H2A%uooo + Q
                                deallocate(Q2)
                                allocate(Q2(nua,noa,noa,noa))
                                call einsum('amie,ej->amij',vA%uoou,t1a,Q)
                                H2A%uooo = H2A%uooo + Q
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q2)
                                H2A%uooo = H2A%uooo - Q2
                                deallocate(Q2)
                                allocate(Q2(noa,noa,nua,noa))
                                call einsum('nmef,fj->nmej',vA%oouu,t1a,Q2)
                                call einsum('nmej,an->amej',Q2,t1a,Q3)
                                call einsum('amej,ei->amij',Q3,t1a,Q)
                                H2A%uooo = H2A%uooo - 0.5*Q
                                deallocate(Q,Q2,Q3)

                                ! H2A(uuou)
                                H2A%uuou = vA%uuou
                                allocate(Q(nua,nua,noa,nua),Q2(nua,noa,noa,nua))
                                call einsum('abfe,fi->abie',vA%uuuu,t1a,Q)
                                H2A%uuou = H2A%uuou + 0.5*Q
                                call einsum('mnie,am->anie',vA%ooou,t1a,Q2)
                                call einsum('anie,bn->abie',Q2,t1a,Q)
                                H2A%uuou = H2A%uuou + Q
                                deallocate(Q,Q2)

                                ! H2A(uoou)
                                H2A%uoou = vA%uoou
                                allocate(Q(nua,noa,noa,nua),Q2(noa,noa,noa,nua))
                                call einsum('nmfe,fi->nmie',vA%oouu,t1a,Q2)
                                call einsum('nmie,an->amie',vA%ooou+Q2,t1a,Q)
                                H2A%uoou = H2A%uoou - Q
                                call einsum('amfe,fi->amie',vA%uouu,t1a,Q)
                                H2A%uoou = H2A%uoou + Q
                                deallocate(Q,Q2)

                                ! H2B(oooo)
                                H2B%oooo = vB%oooo
                                allocate(Q(noa,nob,noa,nob),Q2(noa,nob,noa,nub))
                                call einsum('mnej,ei->mnij',vB%oouo,t1a,Q)
                                H2B%oooo = H2B%oooo + Q
                                call einsum('mnef,ei->mnif',vB%oouu,t1a,Q2)
                                call einsum('mnif,fj->mnij',vB%ooou+Q2,t1b,Q)
                                H2B%oooo = H2B%oooo + Q
                                deallocate(Q,Q2)

                                ! H2B(vvvv)
                                H2B%uuuu = vB%uuuu
                                allocate(Q(nua,nub,nua,nub),Q2(nua,nob,nua,nub))
                                call einsum('mbef,am->abef',vB%ouuu,t1a,Q)
                                H2B%uuuu = H2B%uuuu - Q
                                call einsum('mnef,am->anef',vB%oouu,t1a,Q2)
                                call einsum('anef,bn->abef',vB%uouu-Q2,t1b,Q)
                                H2B%uuuu = H2B%uuuu - Q
                                deallocate(Q,Q2)

                                ! H2B(uoou)
                                H2B%uoou = vB%uoou
                                allocate(Q(nua,nob,noa,nub),Q2(noa,nob,noa,nub))
                                call einsum('nmfe,fi->nmie',vB%oouu,t1a,Q2)
                                call einsum('nmie,an->amie',vB%ooou+Q2,t1a,Q)
                                H2B%uoou = H2B%uoou - Q
                                call einsum('amfe,fi->amie',vB%uouu,t1a,Q)
                                H2B%uoou = H2B%uoou + Q
                                deallocate(Q,Q2)
        
                                ! H2B(ouou)
                                H2B%ouou = vB%ouou 
                                allocate(Q(noa,nub,noa,nub),Q2(noa,nub,nua,nub))
                                call einsum('mnfe,an->mafe',vB%oouu,t1b,Q2)
                                call einsum('mafe,fi->maie',vB%ouuu-Q2,t1a,Q)
                                H2B%ouou = H2B%ouou + Q
                                call einsum('mnie,an->maie',vB%ooou,t1b,Q)
                                H2B%ouou = H2B%ouou - Q
                                deallocate(Q,Q2)

                                ! H2B(uouo)
                                H2B%uouo = vB%uouo
                                allocate(Q(nua,nob,nua,nob),Q2(noa,nob,nua,nob))
                                call einsum('nmef,fi->nmei',vB%oouu,t1b,Q2)
                                call einsum('nmei,an->amei',vB%oouo+Q2,t1a,Q)
                                H2B%uouo = H2B%uouo - Q
                                call einsum('amef,fi->amei',vB%uouu,t1b,Q)    
                                H2B%uouo = H2B%uouo + Q
                                deallocate(Q,Q2)
          
                               ! H2B(ouuo)
                               H2B%ouuo = vB%ouuo
                               allocate(Q(noa,nub,nua,nob),Q2(noa,nob,nua,nob))
                               call einsum('maef,fi->maei',vB%ouuu,t1b,Q)
                               H2B%ouuo = H2B%ouuo + Q
                               call einsum('mnef,fi->mnei',vB%oouu,t1b,Q2)
                               call einsum('mnei,an->maei',vB%oouo+Q2,t1b,Q)
                               H2B%ouuo = H2B%ouuo - Q
                               deallocate(Q,Q2)

                               ! H2B(ouoo)
                               H2B%ouoo = vB%ouoo
                               allocate(Q(noa,nub,noa,nob),Q2(noa,nob,noa,nob),Q3(noa,nob,noa,nob))
                               call einsum('mnif,fj->mnij',vB%ooou,t1b,Q2)
                               Q2 = Q2 + vB%oooo
                               call einsum('mnfj,fi->mnij',vB%oouo,t1a,Q3)
                               Q2 = Q2 + Q3
                               call einsum('mnij,bn->mbij',Q2,t1b,Q)
                               H2B%ouoo = H2B%ouoo - Q
                               deallocate(Q2,Q3)
                               allocate(Q2(noa,nub,nua,nob))
                               call einsum('mbef,fj->mbej',vB%ouuu,t1b,Q2)
                               Q2 = Q2 + vB%ouuo
                               call einsum('mbej,ei->mbij',Q2,t1a,Q)
                               H2B%ouoo = H2B%ouoo + Q
                               deallocate(Q,Q2)
    
                              ! H2B(uooo)
                              H2B%uooo = vB%uooo
                              allocate(Q(nua,nob,noa,nob),Q2(nua,nob,noa,nub),Q3(nua,nob,nua,nub))
                              call einsum('nmfe,an->amfe',vB%oouu,t1a,Q3)
                              Q3 = -1.0*Q3 + vB%uouu
                              call einsum('amfe,fi->amie',Q3,t1a,Q2)
                              Q2 = Q2 + vB%uoou
                              call einsum('amie,ej->amij',Q2,t1b,Q)
                              H2B%uooo = vB%uooo + Q
                              deallocate(Q,Q2,Q3)
 
                              ! H2B(uuuo)
                              H2B%uuuo = vB%uuuo
                              allocate(Q(nua,nub,nua,nob))
                              call einsum('abef,fj->abej',vB%uuuu,t1b,Q)
                              H2B%uuuo = H2B%uuuo + Q
                              call einsum('anej,bn->abej',vB%uouo,t1b,Q)
                              H2B%uuuo = H2B%uuuo - Q
                              deallocate(Q)
                         
                              ! H2B(uuou)
                              H2B%uuou = vB%uuou
                              allocate(Q(nua,nub,noa,nub))
                              call einsum('mbie,am->abie',vB%ouou,t1a,Q)
                              H2B%uuou = H2B%uuou - Q
                              deallocate(Q)
 
                              ! H2C(oooo)
                              H2C%oooo = vC%oooo
                              allocate(Q(nob,nob,nob,nob),Q2(nob,nob,nub,nob),Q3(nob,nob,nob,nob))
                              call einsum('mnef,fj->mnej',vC%oouu,t1b,Q2)
                              call einsum('mnej,ei->mnij',vC%oouo+0.5*Q2,t1b,Q)
                              H2C%oooo = H2C%oooo + Q
                              call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q3)
                              H2C%oooo = H2C%oooo - Q3
                              deallocate(Q,Q2,Q3)
                               
                              ! H2C(uuuu)
                              H2C%uuuu = vC%uuuu
                              allocate(Q(nub,nub,nub,nub),Q2(nob,nub,nub,nub),Q3(nub,nub,nub,nub))
                              call einsum('mnef,bn->mbef',vC%oouu,t1b,Q2)
                              call einsum('mbef,am->abef',vC%ouuu-0.5*Q2,t1b,Q)
                              H2C%uuuu = H2C%uuuu - Q
                              call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q3)
                              H2C%uuuu = H2C%uuuu + Q3
                              deallocate(Q,Q2,Q3)
 
                              ! H2C(uoou)
                              H2C%uoou = vC%uoou
                              allocate(Q(nub,nob,nob,nub),Q2(nob,nob,nob,nub))
                              call einsum('nmfe,fi->nmie',vC%oouu,t1b,Q2)
                              call einsum('nmie,an->amie',vC%ooou+Q2,t1b,Q)
                              H2C%uoou = H2C%uoou - Q
                              call einsum('amfe,fi->amie',vC%uouu,t1b,Q)
                              H2C%uoou = H2C%uoou + Q
                              deallocate(Q,Q2)
 
                              ! H2C(uooo)
                              H2C%uooo = vC%uooo
                              allocate(Q(nub,nob,nob,nob),Q2(nub,nob,nob,nub),Q3(nub,nob,nub,nob))
                              call einsum('nmij,an->amij',vC%oooo,t1b,Q)
                              H2C%uooo = H2C%uooo - 0.5*Q
                              call einsum('amef,ei->amif',vC%uouu,t1b,Q2)
                              call einsum('amif,fj->amij',Q2,t1b,Q)
                              H2C%uooo = H2C%uooo + Q
                              deallocate(Q2)
                              allocate(Q2(nub,nob,nob,nob))
                              call einsum('amie,ej->amij',vC%uoou,t1b,Q)
                              H2C%uooo = H2C%uooo + Q
                              call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q2)
                              H2C%uooo = H2C%uooo - Q2
                              deallocate(Q2)
                              allocate(Q2(nob,nob,nub,nob))
                              call einsum('nmef,fj->nmej',vC%oouu,t1b,Q2)
                              call einsum('nmej,an->amej',Q2,t1b,Q3)
                              call einsum('amej,ei->amij',Q3,t1b,Q)
                              H2C%uooo = H2C%uooo - 0.5*Q
                              deallocate(Q,Q2,Q3)
                              
 
                              ! H2C(uuuo)
                              H2C%uuuo = vC%uuuo
                              allocate(Q(nub,nub,nub,nob),Q2(nub,nob,nub,nob))
                              call einsum('abef,fj->abej',vC%uuuu,t1b,Q)
                              H2C%uuuo = H2C%uuuo + 0.5*Q
                              call einsum('mnej,am->anej',vC%oouo,t1b,Q2)
                              call einsum('anej,bn->abej',Q2,t1b,Q)
                              H2C%uuuo = H2C%uuuo + Q
                              deallocate(Q,Q2)
 
                        end associate


                end subroutine hbar_ccs_intermediates

end module ccsd_module

