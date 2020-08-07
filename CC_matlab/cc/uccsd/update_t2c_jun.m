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
    q10 = einsum_kg(sys.vA_oovv,t2c,'mnef,efjn->mj');
    q15 = einsum_kg(sys.vA_oovv,t1b,'mnef,fn->em');
    q16 = einsum_kg(q15,t1b,'em,bm->eb');
    s25 = einsum_kg(sys.vA_ooov,t2c,'mnie,aejm->inja');

    x1 = np.einsum('ijmb->ijbm', v_a[o,o,o,u]) - np.einsum('ibmj->ijbm', s49) + 0.5*np.einsum('bmij->ijbm', s30)
    z1 = np.einsum('ijbm,am->ijba', x1, t1b, optimize=True)
    z2 = np.einsum('ieab,ej->iabj', v_a[o,u,u,u], t1b, optimize=True)
    x2 = np.einsum('mbei->ibem', v_b[o,u,u,o]) + np.einsum('bemi->ibem', s19) + 0.5*np.einsum('emib->ibem', s33)
    z3 = np.einsum('ibem,eamj->ibja', x2, t2b, optimize=True)
    x3 = np.einsum('mi->im', f[o,o]) + np.einsum('im->im', q1) + np.einsum('mi->im', q3) + np.einsum('im->im', q5) + np.einsum('mi->im', q12) + np.einsum('mi->im', q7) + np.einsum('im->im', q14) + 0.5*np.einsum('mi->im', q10)
    z4 = np.einsum('im,abjm->ijab', x3, t2c, optimize=True)
    x4 = np.einsum('be->be', f[u,u]) + np.einsum('be->be', q2) - np.einsum('eb->be', q4) - np.einsum('be->be', q6) - np.einsum('eb->be', q8) + 0.5*np.einsum('eb->be', q9) - np.einsum('eb->be', q13) - np.einsum('eb->be', q16)
    z5 = np.einsum('be,aeij->bija', x4, t2c, optimize=True)
    x5 = np.einsum('mnij->ijmn', v_a[o,o,o,o]) + 0.5*np.einsum('mnij->ijmn', s40) + np.einsum('imnj->ijmn', s54)
    z6 = np.einsum('ijmn,abmn->ijab', x5, t2c, optimize=True)
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
