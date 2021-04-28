module ccsd_module

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use einsum_module, only: einsum
        use permutils, only: reorder_stripe

        implicit none

        contains

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

                end subroutine update_t1b


                subroutine update_t2a(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(in) :: shift
                        real, intent(inout) :: t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a)

                end subroutine update_t2a


                subroutine update_t2b(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(in) :: shift
                        real, intent(inout) :: t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b)
               
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

