module hbar

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use einsum_module, only: einsum
        use permutils, only: reorder_stripe

        implicit none

        contains

                subroutine hbar_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,H1A,H1B,H2A,H2B,H2C)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        type(e1int_t), intent(out) :: H1A, H1B
                        type(e2int_t), intent(out) :: H2A, H2B, H2C
                        real :: I2A_vovv(sys%Nunocc_a,sys%Nocc_a,sys%Nunocc_a,sys%Nunocc_a), &
                                I2A_ooov(sys%Nocc_a,sys%Nocc_a,sys%Nocc_a,sys%Nunocc_a), &
                                I2B_vovv(sys%Nunocc_a,sys%Nocc_b,sys%Nunocc_a,sys%Nunocc_b), &
                                I2B_ooov(sys%Nocc_a,sys%Nocc_b,sys%Nocc_a,sys%Nunocc_b), &
                                I2B_ovvv(sys%Nocc_a,sys%Nunocc_b,sys%Nunocc_a,sys%Nunocc_b), &
                                I2B_oovo(sys%Nocc_a,sys%Nocc_b,sys%Nunocc_a,sys%Nocc_b), &
                                I2C_vovv(sys%Nunocc_b,sys%Nocc_b,sys%Nunocc_b,sys%Nunocc_b), &
                                I2C_ooov(sys%Nocc_b,sys%Nocc_b,sys%Nocc_b,sys%Nunocc_b)
                        real, allocatable :: Z(:,:), Q(:,:,:,:), Q1(:,:,:,:), Q2(:,:,:,:)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

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
                                call einsum('je,ei->ji',H1A%ou,t1a,Z)
                                H1A%oo = H1A%oo + Z
                                call einsum('jmie,em->ji',vA%ooou,t1a,Z)
                                H1A%oo = H1A%oo + Z
                                call einsum('jmie,em->ji',vB%ooou,t1b,Z)
                                H1A%oo = H1A%oo + Z
                                call einsum('jnef,efin->ji',vA%oouu,t2a,Z)
                                H1A%oo = H1A%oo + 0.5*Z
                                call einsum('jnef,efin->ji',vB%oouu,t2b,Z)
                                H1A%oo = H1A%oo + Z
                                deallocate(Z)

                                ! H1B(oo)
                                H1B%oo = fB%oo
                                allocate(Z(nob,nob))
                                call einsum('je,ei->ji',H1B%ou,t1b,Z)
                                H1B%oo = H1B%oo + Z
                                call einsum('jmie,em->ji',vC%ooou,t1b,Z)
                                H1B%oo = H1B%oo + Z
                                call einsum('mjei,em->ji',vB%oouo,t1a,Z)
                                H1B%oo = H1B%oo + Z
                                call einsum('jnef,efin->ji',vC%oouu,t2c,Z)
                                H1B%oo = H1B%oo + 0.5*Z
                                call einsum('njfe,feni->ji',vB%oouu,t2b,Z)
                                H1B%oo = H1B%oo + Z
                                deallocate(Z)

                                ! H1A(vv)
                                H1A%uu = fA%uu
                                allocate(Z(nua,nua))
                                call einsum('mb,am->ab',H1A%ou,t1a,Z)
                                H1A%uu = H1A%uu - Z
                                call einsum('ambe,em->ab',vA%uouu,t1a,Z)
                                H1A%uu = H1A%uu + Z
                                call einsum('ambe,em->ab',vB%uouu,t1b,Z)
                                H1A%uu = H1A%uu + Z
                                call einsum('mnbf,afmn->ab',vA%oouu,t2a,Z)
                                H1A%uu = H1A%uu - 0.5*Z
                                call einsum('mnbf,afmn->ab',vB%oouu,t2b,Z)
                                H1A%uu = H1A%uu - Z
                                deallocate(Z)

                                ! H1B(vv)
                                H1B%uu = fB%uu
                                allocate(Z(nub,nub))
                                call einsum('mb,am->ab',H1B%ou,t1b,Z)
                                H1B%uu = H1B%uu - Z
                                call einsum('ambe,em->ab',vC%uouu,t1b,Z)
                                H1B%uu = H1B%uu + Z
                                call einsum('maeb,em->ab',vB%ouuu,t1a,Z)
                                H1B%uu = H1B%uu + Z
                                call einsum('mnbf,afmn->ab',vC%oouu,t2c,Z)
                                H1B%uu = H1B%uu - 0.5*Z
                                call einsum('nmfb,fanm->ab',vB%oouu,t2b,Z)
                                H1B%uu = H1B%uu - Z
                                deallocate(Z)

                                ! H2A(uouu)
                                allocate(Q(nua,noa,nua,nua))
                                call einsum('mnfe,an->amef',vA%oouu,t1a,Q)
                                I2A_vovv = vA%uouu - 0.5*Q
                                H2A%uouu = I2A_vovv - 0.5*Q
                                deallocate(Q)

                                ! H2A(ooou)
                                allocate(Q(noa,noa,noa,nua))
                                call einsum('mnfe,fi->mnie',vA%oouu,t1a,Q)
                                I2A_ooov = vA%ooou + 0.5*Q
                                H2A%ooou = I2A_ooov + 0.5*Q
                                deallocate(Q)

                                ! H2B(uouu)
                                allocate(Q(nua,nob,nua,nub))
                                call einsum('nmef,an->amef',vB%oouu,t1a,Q)
                                I2B_vovv = vB%uouu - 0.5*Q
                                H2B%uouu = I2B_vovv - 0.5*Q
                                deallocate(Q)

                                ! H2B(ooou)
                                allocate(Q(noa,nob,noa,nub))
                                call einsum('mnfe,fi->mnie',vB%oouu,t1a,Q)
                                I2B_ooov = vB%ooou + 0.5*Q
                                H2B%ooou = I2B_ooov + 0.5*Q
                                deallocate(Q)

                                ! H2B(ouuu)
                                allocate(Q(noa,nub,nua,nub))
                                call einsum('mnef,an->maef',vB%oouu,t1b,Q)
                                I2B_ovvv = vB%ouuu - 0.5*Q
                                H2B%ouuu = I2B_ovvv - 0.5*Q
                                deallocate(Q)

                                ! H2B(oouo)
                                allocate(Q(noa,nob,nua,nob))
                                call einsum('nmef,fi->nmei',vB%oouu,t1b,Q)
                                I2B_oovo = vB%oouo + 0.5*Q
                                H2B%oouo = I2B_oovo + 0.5*Q
                                deallocate(Q)

                                ! H2C(uouu)
                                allocate(Q(nub,nob,nub,nub))
                                call einsum('mnfe,an->amef',vC%oouu,t1b,Q)
                                I2C_vovv = vC%uouu - 0.5*Q
                                H2C%uouu = I2C_vovv - 0.5*Q
                                deallocate(Q)

                                ! H2C(ooou)
                                allocate(Q(nob,nob,nob,nub))
                                call einsum('mnfe,fi->mnie',vC%oouu,t1b,Q)
                                I2C_ooov = vC%ooou + 0.5*Q
                                H2C%ooou = I2C_ooov + 0.5*Q
                                deallocate(Q)

                                ! H2A(uuuu)
                                H2A%uuuu = vA%uuuu
                                allocate(Q(nua,nua,nua,nua),Q1(nua,nua,nua,nua))
                                call einsum('bmfe,am->abef',I2A_vovv,t1a,Q)
                                Q = -1.0*Q
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                Q = Q - Q1
                                H2A%uuuu = H2A%uuuu + Q
                                call einsum('mnef,abmn->abef',vA%oouu,t2a,Q)
                                H2A%uuuu = H2A%uuuu + 0.5*Q
                                deallocate(Q,Q1)

                                ! H2B(uuuu)
                                H2B%uuuu = vB%uuuu
                                allocate(Q(nua,nub,nua,nub))
                                call einsum('mbef,am->abef',I2B_ovvv,t1a,Q)
                                H2B%uuuu = H2B%uuuu - Q
                                call einsum('amef,bm->abef',I2B_vovv,t1b,Q)
                                H2B%uuuu = H2B%uuuu - Q
                                call einsum('mnef,abmn->abef',vB%oouu,t2b,Q)
                                H2B%uuuu = H2B%uuuu + Q
                                deallocate(Q)

                                ! H2C(uuuu) 
                                H2C%uuuu = vC%uuuu
                                allocate(Q(nub,nub,nub,nub),Q1(nub,nub,nub,nub))
                                call einsum('bmfe,am->abef',I2C_vovv,t1b,Q)
                                Q = -1.0*Q
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                Q = Q - Q1
                                H2C%uuuu = H2C%uuuu + Q
                                call einsum('mnef,abmn->abef',vC%oouu,t2c,Q)
                                H2C%uuuu = H2C%uuuu + 0.5*Q
                                deallocate(Q,Q1)

                                ! H2A(oooo)
                                H2A%oooo = vA%oooo
                                allocate(Q(noa,noa,noa,noa),Q1(noa,noa,noa,noa))
                                call einsum('nmje,ei->mnij',I2A_ooov,t1a,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                Q = Q - Q1
                                H2A%oooo = H2A%oooo + Q
                                call einsum('mnef,efij->mnij',vA%oouu,t2a,Q)
                                H2A%oooo = H2A%oooo + 0.5*Q
                                deallocate(Q,Q1)

                                ! H2B(oooo)
                                H2B%oooo = vB%oooo
                                allocate(Q(noa,nob,noa,nob))
                                call einsum('mnej,ei->mnij',I2B_oovo,t1a,Q)
                                H2B%oooo = H2B%oooo + Q
                                call einsum('mnie,ej->mnij',I2B_ooov,t1b,Q)
                                H2B%oooo = H2B%oooo + Q
                                call einsum('mnef,efij->mnij',vB%oouu,t2b,Q)
                                H2B%oooo = H2B%oooo + Q
                                deallocate(Q)

                                ! H2C(oooo)
                                H2C%oooo = vC%oooo
                                allocate(Q(nob,nob,nob,nob),Q1(nob,nob,nob,nob))
                                call einsum('nmje,ei->mnij',I2C_ooov,t1b,Q)
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                Q = Q - Q1
                                H2C%oooo = H2C%oooo + Q
                                call einsum('mnef,efij->mnij',vC%oouu,t2c,Q)
                                H2C%oooo = H2C%oooo + 0.5*Q
                                deallocate(Q,Q1)

                                ! H2A(uoou)
                                H2A%uoou = vA%uoou
                                allocate(Q(nua,noa,noa,nua))
                                call einsum('amfe,fi->amie',I2A_vovv,t1a,Q)
                                H2A%uoou = H2A%uoou + Q
                                call einsum('nmie,an->amie',I2A_ooov,t1a,Q)
                                H2A%uoou = H2A%uoou - Q
                                call einsum('nmfe,afin->amie',vA%oouu,t2a,Q)
                                H2A%uoou = H2A%uoou + Q
                                call einsum('mnef,afin->amie',vB%oouu,t2b,Q)
                                H2A%uoou = H2A%uoou + Q
                                deallocate(Q)

                                ! H2B(uoou)
                                H2B%uoou = vB%uoou
                                allocate(Q(nua,noa,noa,nua))
                                call einsum('amfe,fi->amie',I2B_vovv,t1a,Q)
                                H2B%uoou = H2B%uoou + Q
                                call einsum('nmie,an->amie',I2B_ooov,t1a,Q)
                                H2B%uoou = H2B%uoou - Q
                                call einsum('nmfe,afin->amie',vB%oouu,t2a,Q)
                                H2B%uoou = H2B%uoou + Q
                                call einsum('nmfe,afin->amie',vC%oouu,t2b,Q)
                                H2B%uoou = H2B%uoou + Q
                                deallocate(Q)

                                ! H2B(ouuo)
                                H2B%ouuo = vB%ouuo
                                allocate(Q(noa,nub,nua,nob))
                                call einsum('maef,fi->maei',I2B_ovvv,t1b,Q)
                                H2B%ouuo = H2B%ouuo + Q
                                call einsum('mnei,an->maei',I2B_oovo,t1b,Q)
                                H2B%ouuo = H2B%ouuo - Q
                                call einsum('mnef,afin->maei',vB%oouu,t2c,Q)
                                H2B%ouuo = H2B%ouuo + Q
                                call einsum('mnef,fani->maei',vA%oouu,t2b,Q)
                                H2B%ouuo = H2B%ouuo + Q
                                deallocate(Q)

                                ! H2B(ouou)
                                H2B%ouou = vB%ouou
                                allocate(Q(noa,nub,noa,nub))
                                call einsum('mafe,fi->maie',I2B_ovvv,t1a,Q)
                                H2B%ouou = H2B%ouou + Q
                                call einsum('mnie,an->maie',I2B_ooov,t1b,Q)
                                H2B%ouou = H2B%ouou - Q
                                call einsum('mnfe,fain->maie',vB%oouu,t2b,Q)
                                H2B%ouou = H2B%ouou - Q
                                deallocate(Q)

                                ! H2B(uouo)
                                H2B%uouo = vB%uouo
                                allocate(Q(nua,nob,nua,nob))
                                call einsum('nmei,an->amei',I2B_oovo,t1a,Q)
                                H2B%uouo = H2B%uouo - Q
                                call einsum('amef,fi->amei',I2B_vovv,t1b,Q)
                                H2B%uouo = H2B%uouo + Q
                                call einsum('nmef,afni->amei',vB%oouu,t2b,Q)
                                H2B%uouo = H2B%uouo - Q
                                deallocate(Q)

                                ! H2C(uoou)
                                H2C%uoou = vC%uoou
                                allocate(Q(nub,nob,nob,nub))
                                call einsum('amfe,fi->amie',I2C_vovv,t1b,Q)
                                H2C%uoou = H2C%uoou + Q
                                call einsum('nmie,an->amie',I2C_ooov,t1b,Q)
                                H2C%uoou = H2C%uoou - Q
                                call einsum('nmfe,afin->amie',vC%oouu,t2c,Q)
                                H2C%uoou = H2C%uoou + Q
                                call einsum('nmfe,fani->amie',vB%oouu,t2b,Q)
                                H2C%uoou = H2C%uoou + Q
                                deallocate(Q)

                                ! H2A(uooo)
                                H2A%uooo = vA%uooo
                                allocate(Q(nua,noa,noa,noa),Q1(nua,noa,noa,noa),Q2(nua,noa,noa,nua))
                                call einsum('mnjf,afin->amij',H2A%ooou,t2a,Q)
                                call einsum('mnjf,afin->amij',H2B%ooou,t2b,Q1)
                                Q = Q + Q1
                                call einsum('amef,ei->amif',vA%uouu,t1a,Q2)
                                Q2 = vA%uoou + 0.5*Q2
                                call einsum('amif,fj->amij',Q2,t1a,Q1)
                                Q = Q + Q1
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                Q = Q - Q1
                                H2A%uooo = H2A%uooo + Q
                                call einsum('me,aeij->amij',H1A%ou,t2a,Q)
                                H2A%uooo = H2A%uooo + Q
                                call einsum('nmij,an->amij',H2A%oooo,t1a,Q)
                                H2A%uooo = H2A%uooo - Q
                                call einsum('amef,efij->amij',vA%uouu,t2a,Q)
                                H2A%uooo = H2A%uooo + 0.5*Q
                                deallocate(Q,Q1,Q2)

                                ! H2B(uooo)
                                H2B%uooo = vB%uooo
                                allocate(Q(nua,nob,noa,nob),Q1(nua,nob,noa,nub))
                                call einsum('amfe,fi->amie',vB%uouu,t1a,Q1)
                                Q1 = vB%uoou + Q1
                                call einsum('amie,ej->amij',Q1,t1b,Q)
                                H2B%uooo = H2B%uooo + Q
                                call einsum('me,aeij->amij',H1B%ou,t2b,Q)
                                H2B%uooo = H2B%uooo + Q
                                call einsum('nmij,an->amij',H2B%oooo,t1a,Q)
                                H2B%uooo = H2B%uooo - Q
                                call einsum('mnjf,afin->amij',H2C%ooou,t2b,Q)
                                H2B%uooo = H2B%uooo + Q
                                call einsum('nmfj,afin->amij',H2B%oouo,t2a,Q)
                                H2B%uooo = H2B%uooo + Q
                                call einsum('nmif,afnj->amij',H2B%ooou,t2b,Q)
                                H2B%uooo = H2B%uooo - Q
                                call einsum('amej,ei->amij',vB%uouo,t1a,Q)
                                H2B%uooo = H2B%uooo + Q
                                call einsum('amef,efij->amij',vB%uouu,t2b,Q)
                                H2B%uooo = H2B%uooo + Q
                                deallocate(Q,Q1)

                                ! H2B(ouoo)
                                H2B%ouoo = vB%ouoo
                                allocate(Q(noa,nub,noa,nob),Q1(noa,nub,noa,nub))
                                call einsum('mafe,fj->maje',vB%ouuu,t1a,Q1)
                                Q1 = vB%ouou + Q1
                                call einsum('maje,ei->maji',Q1,t1b,Q)
                                H2B%ouoo = H2B%ouoo + Q
                                call einsum('me,eaji->maji',H1A%ou,t2b,Q)
                                H2B%ouoo = H2B%ouoo + Q
                                call einsum('mnji,an->maji',H2B%oooo,t1b,Q)
                                H2B%ouoo = H2B%ouoo - Q
                                call einsum('mnjf,fani->maji',H2A%ooou,t2b,Q)
                                H2B%ouoo = H2B%ouoo + Q
                                call einsum('mnjf,fani->maji',H2B%ooou,t2c,Q)
                                H2B%ouoo = H2B%ouoo + Q
                                call einsum('mnfi,fajn->maji',H2B%oouo,t2b,Q)
                                H2B%ouoo = H2B%ouoo - Q
                                call einsum('maei,ej->maji',vB%ouuo,t1a,Q)
                                H2B%ouoo = H2B%ouoo + Q
                                call einsum('mafe,feji->maji',vB%ouuu,t2b,Q)
                                H2B%ouoo = H2B%ouoo + Q
                                deallocate(Q,Q1)

                                ! H2C(uooo)
                                H2C%uooo = vC%uooo
                                allocate(Q(nub,nob,nob,nob),Q1(nub,nob,nob,nob),Q2(nub,nob,nob,nub))
                                call einsum('mnjf,afin->amij',H2C%ooou,t2c,Q)
                                call einsum('nmfj,fani->amij',H2B%oouo,t2b,Q1)
                                Q = Q + Q1
                                call einsum('amef,ei->amif',vC%uouu,t1b,Q2)
                                Q2 = vC%uoou + 0.5*Q2
                                call einsum('amif,fj->amij',Q2,t1b,Q1)
                                Q = Q + Q1
                                call reorder_stripe(4,shape(Q),size(Q),'1243',Q,Q1)
                                Q = Q - Q1
                                H2C%uooo = H2C%uooo + Q
                                call einsum('me,aeij->amij',H1B%ou,t2c,Q)
                                H2C%uooo = H2C%uooo + Q
                                call einsum('nmij,an->amij',H2C%oooo,t1b,Q)
                                H2C%uooo = H2C%uooo - Q
                                call einsum('amef,efij->amij',vC%uouu,t2c,Q)
                                H2C%uooo = H2C%uooo + 0.5*Q
                                deallocate(Q,Q1,Q2)

                                ! H2A(uuou)
                                H2A%uuou = vA%uuou
                                allocate(Q(nua,nua,noa,nua),Q1(nua,nua,noa,nua),Q2(noa,nua,noa,nua))
                                call einsum('bnef,afin->abie',H2A%uouu,t2a,Q)
                                call einsum('bnef,afin->abie',H2B%uouu,t2b,Q1)
                                Q = Q + Q1
                                call einsum('mnie,bn->mbie',vA%ooou,t1a,Q2)
                                Q2 = vA%ouou - 0.5*Q2
                                call einsum('mbie,am->abie',Q2,t1a,Q1)
                                Q = Q - Q1
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                Q = Q - Q1
                                H2A%uuou = H2A%uuou + Q
                                call einsum('me,abim->abie',H1A%ou,t2a,Q)
                                H2A%uuou = H2A%uuou - Q
                                call einsum('abfe,fi->abie',H2A%uuuu,t1a,Q)
                                H2A%uuou = H2A%uuou + Q
                                call einsum('mnie,abmn->abie',vA%ooou,t2a,Q)
                                H2A%uuou = H2A%uuou + 0.5*Q
                                deallocate(Q,Q1,Q2)

                                ! H2B(uuou)
                                H2B%uuou = vB%uuou
                                allocate(Q(nua,nub,noa,nub),Q1(noa,nub,noa,nub))
                                call einsum('mnie,bn->mbie',vB%ooou,t1b,Q1)
                                Q1 = vB%ouou - Q1
                                call einsum('mbie,am->abie',Q1,t1a,Q)
                                H2B%uuou = H2B%uuou - Q
                                call einsum('me,abim->abie',H1B%ou,t2b,Q)
                                H2B%uuou = H2B%uuou - Q
                                call einsum('abfe,fi->abie',H2B%uuuu,t1a,Q)
                                H2B%uuou = H2B%uuou + Q
                                call einsum('nbfe,afin->abie',H2B%ouuu,t2a,Q)
                                H2B%uuou = H2B%uuou + Q
                                call einsum('bnef,afin->abie',H2C%uouu,t2b,Q)
                                H2B%uuou = H2B%uuou + Q
                                call einsum('amfe,fbim->abie',H2B%uouu,t2b,Q)
                                H2B%uuou = H2B%uuou - Q
                                call einsum('amie,bm->abie',vB%uoou,t1b,Q)
                                H2B%uuou = H2B%uuou - Q
                                call einsum('nmie,abnm->abie',vB%ooou,t2b,Q)
                                H2B%uuou = H2B%uuou + Q
                                deallocate(Q,Q1)

                                ! H2B(uuuo)
                                H2B%uuuo = vB%uuuo
                                allocate(Q(nua,nub,nua,nob),Q1(nua,nob,nua,nob))
                                call einsum('nmei,bn->bmei',vB%oouo,t1a,Q1)
                                Q1 = vB%uouo - Q1
                                call einsum('bmei,am->baei',Q1,t1b,Q)
                                H2B%uuuo = H2B%uuuo - Q
                                call einsum('me,bami->baei',H1A%ou,t2b,Q)
                                H2B%uuuo = H2B%uuuo - Q
                                call einsum('baef,fi->baei',H2B%uuuu,t1b,Q)
                                H2B%uuuo = H2B%uuuo + Q
                                call einsum('bnef,fani->baei',H2A%uouu,t2b,Q)
                                H2B%uuuo = H2B%uuuo + Q
                                call einsum('bnef,fani->baei',H2B%uouu,t2c,Q)
                                H2B%uuuo = H2B%uuuo + Q
                                call einsum('maef,bfmi->baei',H2B%ouuu,t2b,Q)
                                H2B%uuuo = H2B%uuuo - Q
                                call einsum('naei,bn->baei',vB%ouuo,t1a,Q)
                                H2B%uuuo = H2B%uuuo - Q
                                call einsum('nmei,banm->baei',vB%oouo,t2b,Q)
                                H2B%uuuo = H2B%uuuo + Q
                                deallocate(Q,Q1)

                                ! H2C(uuou)
                                H2C%uuou = vC%uuou
                                allocate(Q(nub,nub,nob,nub),Q1(nub,nub,nob,nub),Q2(nob,nub,nob,nub))
                                call einsum('bnef,afin->abie',H2C%uouu,t2c,Q)
                                call einsum('nbfe,fani->abie',H2B%ouuu,t2b,Q1)
                                Q = Q + Q1
                                call einsum('mnie,bn->mbie',vC%ooou,t1b,Q2)
                                Q2 = vC%ouou - 0.5*Q2
                                call einsum('mbie,am->abie',Q2,t1b,Q1)
                                Q = Q - Q1
                                call reorder_stripe(4,shape(Q),size(Q),'2134',Q,Q1)
                                Q = Q - Q1
                                H2C%uuou = H2C%uuou + Q
                                call einsum('me,abim->abie',H1B%ou,t2c,Q)
                                H2C%uuou = H2C%uuou - Q
                                call einsum('abfe,fi->abie',H2C%uuuu,t1b,Q)
                                H2C%uuou = H2C%uuou + Q
                                call einsum('mnie,abmn->abie',vC%ooou,t2c,Q)
                                H2C%uuou = H2C%uuou + 0.5*Q
                                deallocate(Q,Q1,Q2)

                        end associate


                end subroutine hbar_ccsd

end module hbar
