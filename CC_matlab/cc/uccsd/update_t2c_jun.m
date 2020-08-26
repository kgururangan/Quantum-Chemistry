function [t2c] = update_t2c_jun(t1a,t1b,t2a,t2b,t2c,sys)

    q1 = einsum_kg(sys.vB_oovo,t1a,'mnei,em->in');
    q2 = einsum_kg(sys.vB_ovvv,t1a,'mbef,em->bf');
    s11 = einsum_kg(sys.vA_oooo,t1b,'mnij,am->ijna');
    s13 = einsum_kg(sys.vA_ovov,t1b,'mbie,ej->ibmj');
    s15 = einsum_kg(sys.vA_vvvv,t1b,'abef,ei->abfi');
    s17 = einsum_kg(sys.vB_oovo,t2b,'mnei,eamj->inja');
    s19 = einsum_kg(sys.vB_ovvv,t1b,'mbef,fi->bemi');
    q3 = einsum_kg(sys.fb_ov,t1b,'me,ei->mi');
    q4 = einsum_kg(sys.fb_ov,t1b,'me,bm->eb');
    s23 = einsum_kg(sys.vA_ooov,t1b,'mnie,ej->imnj');
    s47 = einsum_kg(s23,t1b,'imnj,am->ijna');
    q5 = einsum_kg(sys.vA_ooov,t1b,'mnie,en->im');
    s28 = einsum_kg(sys.vA_ovvv,t1b,'mbef,fi->bemi');
    s49 = einsum_kg(s28,t1b,'bfmi,fj->ibmj');
    q6 = einsum_kg(sys.vA_ovvv,t1b,'mbef,fm->be');
    s33 = einsum_kg(sys.vA_oovv,t2b,'mnef,ebmi->fnib');
    s35 = einsum_kg(sys.vB_oovv,t2b,'mnef,ebmi->fnib');
    s51 = einsum_kg(s35,t1b,'fnja,fi->jani');
    q8 = einsum_kg(sys.vB_oovv,t2b,'mnef,ebmn->fb');
    q9 = einsum_kg(sys.vA_oovv,t2c,'mnef,aemn->fa');
    s40 = einsum_kg(sys.vA_oovv,t2c,'mnef,efij->mnij');
    s59 = einsum_kg(s40,t1b,'mnij,am->ijna');
    s43 = einsum_kg(sys.vA_oovv,t2c,'mnef,aejm->fnja');
    s56 = einsum_kg(s43,t1b,'fnja,fi->jani');
    q11 = einsum_kg(sys.vB_oovv,t1a,'mnef,em->fn');
    q12 = einsum_kg(q11,t1b,'fn,fi->ni');
    q13 = einsum_kg(q11,t1b,'fn,bn->fb');
    s30 = einsum_kg(sys.vA_ovvv,t2c,'mbef,efij->bmij');
    q7 = einsum_kg(sys.vB_oovv,t2b,'mnef,efmi->ni');
    s53 = einsum_kg(sys.vA_oovv,t1b,'mnef,ei->fmni');
    q14 = einsum_kg(s53,t1b,'fmni,fn->im');
    s54 = einsum_kg(s53,t1b,'fmni,fj->imnj');
    s62 = einsum_kg(s54,t1b,'imnj,am->ijna');
    q10 = einsum_kg(sys.vC_oovv,t2c,'mnef,efjn->mj');
    q15 = einsum_kg(sys.vA_oovv,t1b,'mnef,fn->em');
    q16 = einsum_kg(q15,t1b,'em,bm->eb');
    s25 = einsum_kg(sys.vA_ooov,t2c,'mnie,aejm->inja');

    x1 = einsum_permute(sys.vC_ooov,'ijmb->ijbm') - einsum_permute(s49,'ibmj->ijbm') + 0.5*einsum_permute(s30,'bmij->ijbm');
    z1 = einsum_kg(x1,t1b,'ijbm,am->ijba');
    z2 = einsum_kg(sys.vA_ovvv,t1b,'ieab,ej->iabj');
    x2 = einsum_permute(sys.vB_ovvo,'mbei->ibem') + einsum_permute('bemi->ibem', s19) + 0.5*einsum_permute('emib->ibem', s33);
    z3 = einsum_kg(x2,t2b,'ibem,eamj->ibja');    
    x3 = einsum_permute('mi->im', sys.fb_oo) + einsum_permute('im->im', q1) + einsum_permute('mi->im', q3)...
         + einsum_permute('im->im', q5) + einsum_permute('mi->im', q12) + einsum_permute('mi->im', q7) ...
         + einsum_permute('im->im', q14) + 0.5*einsum_permute('mi->im', q10);     
    z4 = einsum_kg(x3,t3c,'im,abjm->ijab');
    x4 = einsum_permute('be->be', sys.fb_vv) + einsum_permute('be->be', q2) - einsum_permute('eb->be', q4)...
         - einsum_permute('be->be', q6) - einsum_permute('eb->be', q8) + 0.5*einsum_permute('eb->be', q9)...
         - einsum_permute('eb->be', q13) - einsum_permute('eb->be', q16);     
    z5 = einsum_kg(x4,t2c,'be,aeij->bija');
    x5 = einsum_permute('mnij->ijmn', sys.vC_oooo) + 0.5*einsum_permute('mnij->ijmn', s40) + einsum_permute('imnj->ijmn', s54);
    
    z6 = einsum_kg(x5,t2c,'ijmn,abmn->ijab');
    x6 = np.einsum('mbie->ibem', v_a[o,u,o,u]) - np.einsum('bemi->ibem', s28) - np.einsum('emib->ibem', s35)
    z7 = np.einsum('ibem,aejm->ibja', x6, t2c, optimize=True)
    z8 = np.einsum('abef,efij->abij', v_a[u,u,u,u], t2c, optimize=True)
    x7 = np.einsum('ijna->ijan', s11) + 0.5*np.einsum('ijna->ijan', s59) + np.einsum('ijna->ijan', s62)
    z12 = np.einsum('ijan,bn->ijab', x7, t1b, optimize=True)
    x8 = np.einsum('ibmj->ijbm', s13) + np.einsum('imjb->ijbm', s17) + np.einsum('jbmi->ijbm', s51) + np.einsum('jbmi->ijbm', s56) - np.einsum('imjb->ijbm', s25)
    z14 = np.einsum('ijbm,am->ijba', x8, t1b, optimize=True)
    z16 = np.einsum('abfi,fj->iabj', s15, t1b, optimize=True)
    z24 = np.einsum('imnj,abmn->ijab', s23, t2c, optimize=True)
    z48 = np.einsum('ijna,bn->ijab', s47, t1b, optimize=True)
    z44 = np.einsum('fnja,bfin->jaib', s43, t2c, optimize=True)

    new_t = +np.einsum('abij->abij', v_a[u,u,o,o])
    new_t += -np.einsum('ijba->abij', z1)
    new_t += +np.einsum('ijab->abij', z1)
    new_t += +np.einsum('iabj->abij', z2)
    new_t += -np.einsum('jabi->abij', z2)
    new_t += -np.einsum('ibja->abij', z3)
    new_t += +np.einsum('iajb->abij', z3)
    new_t += +np.einsum('jbia->abij', z3)
    new_t += -np.einsum('jaib->abij', z3)
    new_t += +np.einsum('ijab->abij', z4)
    new_t += -np.einsum('jiab->abij', z4)
    new_t += +np.einsum('bija->abij', z5)
    new_t += -np.einsum('aijb->abij', z5)
    new_t += +0.5*np.einsum('ijab->abij', z6)
    new_t += +np.einsum('ibja->abij', z7)
    new_t += -np.einsum('iajb->abij', z7)
    new_t += -np.einsum('jbia->abij', z7)
    new_t += +np.einsum('jaib->abij', z7)
    new_t += +0.5*np.einsum('abij->abij', z8)
    new_t += +np.einsum('ijab->abij', z12)
    new_t += -np.einsum('ijba->abij', z14)
    new_t += +np.einsum('ijab->abij', z14)
    new_t += +np.einsum('jiba->abij', z14)
    new_t += -np.einsum('jiab->abij', z14)
    new_t += +np.einsum('iabj->abij', z16)
    new_t += +0.5*np.einsum('ijab->abij', z24)
    new_t += -0.5*np.einsum('jiab->abij', z24)
    new_t += -np.einsum('jaib->abij', z44)
    new_t += +np.einsum('iajb->abij', z44)
    new_t += +np.einsum('ijab->abij', z48)
    new_t += -np.einsum('jiab->abij', z48)

    for a = 1:sys.Nvir_beta
        for b = a+1:sys.Nvir_beta
            for i = 1:sys.Nocc_beta
                for j = i+1:sys.Nocc_beta
                    denom = (sys.fb_oo(i,i)+sys.fb_oo(j,j)-sys.fb_vv(a,a)-sys.fb_vv(b,b));
                    coef = new_t(a,b,i,j) / denom;
                    t2c(a,b,i,j) = t2c(a,b,i,j) + coef;               
                    t2c(b,a,i,j) = -t2c(a,b,i,j);
                    t2c(a,b,j,i) = -t2c(a,b,i,j);
                    t2c(b,a,j,i) = t2c(a,b,i,j);
                end
            end
        end
    end

end
