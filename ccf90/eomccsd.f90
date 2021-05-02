module eomccsd

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use permutils, only: reorder_stripe
        use einsum_module, only: einsum
        use cis, only: cis_mat

        implicit none

        contains

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

end module eomccsd
