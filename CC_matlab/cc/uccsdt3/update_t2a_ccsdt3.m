function [t2a] = update_t2a_ccsdt3(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,sys,shift)

    % CCSD part

    s9 = einsum_kg(sys.vA_oooo,t1a,'mnij,am->ijna');
    s11 = einsum_kg(sys.vA_ovov,t1a,'mbie,ej->ibmj');
    s13 = einsum_kg(sys.vA_vvvv,t1a,'abef,ei->abfi');
    q1 = einsum_kg(sys.fa_ov,t1a,'me,ei->mi');
    q2 = einsum_kg(sys.fa_ov,t1a,'me,bm->eb');
    s17 = einsum_kg(sys.vA_ooov,t1a,'mnie,ej->imnj');
    s45 = einsum_kg(s17,t1a,'imnj,am->ijna');
    q3 = einsum_kg(sys.vA_ooov,t1a,'mnie,en->im');
    s22 = einsum_kg(sys.vA_ovvv,t1a,'mbef,fi->bemi');
    s47 = einsum_kg(s22,t1a,'bfmi,fj->ibmj');
    q4 = einsum_kg(sys.vA_ovvv,t1a,'mbef,fm->be');
    s27 = einsum_kg(sys.vB_ooov,t2b,'mnie,aejn->imja');
    s29 = einsum_kg(sys.vB_vovv,t1a,'bmef,ei->bfmi');
    q5 = einsum_kg(sys.vB_ooov,t1b,'mnie,en->im');
    q6 = einsum_kg(sys.vB_vovv,t1b,'bmef,fm->be');
    q7 = einsum_kg(sys.vA_oovv,t2a,'mnef,aemn->fa');
    s34 = einsum_kg(sys.vA_oovv,t2a,'mnef,efij->mnij');
    s55 = einsum_kg(s34,t1a,'mnij,am->ijna');
    s37 = einsum_kg(sys.vA_oovv,t2a,'mnef,aejm->fnja');
    s52 = einsum_kg(s37,t1a,'fnja,fi->jani');
    q10 = einsum_kg(sys.vB_oovv,t2b,'mnef,efjn->mj');
    s41 = einsum_kg(sys.vB_oovv,t2a,'mnef,beim->fnib');
    s43 = einsum_kg(sys.vA_oovv,t2b,'mnef,aejm->fnja');
    s19 = einsum_kg(sys.vA_ooov,t2a,'mnie,aejm->inja');
    s24 = einsum_kg(sys.vA_ovvv,t2a,'mbef,efij->bmij');
    s49 = einsum_kg(sys.vA_oovv,t1a,'mnef,ei->fmni');
    q11 = einsum_kg(s49,t1a,'fmni,fn->im');
    s50 = einsum_kg(s49,t1a,'fmni,fj->imnj');
    s63 = einsum_kg(s50,t1a,'imnj,am->ijna');
    q8 = einsum_kg(sys.vA_oovv,t2a,'mnef,efjn->mj');
    q12 = einsum_kg(sys.vA_oovv,t1a,'mnef,fn->em');
    q13 = einsum_kg(q12,t1a,'em,bm->eb');
    s58 = einsum_kg(sys.vB_oovv,t1a,'mnef,ei->fmni');
    q14 = einsum_kg(s58,t1b,'fmni,fn->im');
    s59 = einsum_kg(s58,t2b,'fmni,afjn->imja');
    q15 = einsum_kg(sys.vB_oovv,t1b,'mnef,fn->em');
    q16 = einsum_kg(q15,t1a,'em,bm->eb');
    q9 = einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ea');

    x1 = einsum_permute(sys.vA_ooov,'ijmb->ijbm') - einsum_permute('ibmj->ijbm', s47) + 0.5*einsum_permute('bmij->ijbm', s24);
    z1 = einsum_kg(x1,t1a,'ijbm,am->ijba');
    z2 = einsum_kg(sys.vA_ovvv,t1a,'ieab,ej->iabj');
    x2 = einsum_permute('mi->im', sys.fa_oo) + einsum_permute('mi->im', q1) + einsum_permute('im->im', q3) + einsum_permute('im->im', q5) + einsum_permute('mi->im', q10) + einsum_permute('im->im', q11) + 0.5*einsum_permute('mi->im', q8) + einsum_permute('im->im', q14);
    z3 = einsum_kg(x2,t2a,'im,abjm->ijab');
    x3 = einsum_permute('be->be', sys.fa_vv) - einsum_permute('eb->be', q2) - einsum_permute('be->be', q4) + einsum_permute('be->be', q6) + 0.5*einsum_permute('eb->be', q7) - einsum_permute('eb->be', q13) - einsum_permute('eb->be', q16) - einsum_permute('eb->be', q9);
    z4 = einsum_kg(x3,t2a,'be,aeij->bija');
    x4 = einsum_permute(sys.vA_oooo,'mnij->ijmn') + 0.5*einsum_permute('mnij->ijmn', s34) + einsum_permute('imnj->ijmn', s50);
    z5 = einsum_kg(x4,t2a,'ijmn,abmn->ijab');
    x5 = einsum_permute(sys.vA_ovov,'mbie->ibem') - einsum_permute('bemi->ibem', s22);
    z6 = einsum_kg(x5,t2a,'ibem,aejm->ibja');
    z7 = einsum_kg(sys.vA_vvvv,t2a,'abef,efij->abij');
    x6 = einsum_permute(sys.vB_ovvo,'iebm->ibem') + einsum_permute('bemi->ibem', s29) + einsum_permute('emib->ibem', s41);
    z8 = einsum_kg(x6,t2b,'ibem,aejm->ibja');
    x7 = einsum_permute('ijna->ijan', s9) + 0.5*einsum_permute('ijna->ijan', s55) + einsum_permute('ijna->ijan', s63);
    z10 = einsum_kg(x7,t1a,'ijan,bn->ijab');
    x8 = einsum_permute('ibmj->ijbm', s11) + einsum_permute('imjb->ijbm', s27) + einsum_permute('jbmi->ijbm', s52) - einsum_permute('imjb->ijbm', s19) + einsum_permute('imjb->ijbm', s59);
    z12 = einsum_kg(x8,t1a,'ijbm,am->ijba');
    z14 = einsum_kg(s13,t1a,'abfi,fj->iabj');
    z18 = einsum_kg(s17,t2a,'imnj,abmn->ijab');
    z46 = einsum_kg(s45,t1a,'ijna,bn->ijab');
    z38 = einsum_kg(s37,t2a,'fnja,bfin->jaib');
    z44 = einsum_kg(s43,t2b,'fnja,bfin->jaib');

    X2A_abij = +einsum_permute('abij->abij', sys.vA_vvoo);
    X2A_abij = X2A_abij -einsum_permute('ijba->abij', z1);
    X2A_abij = X2A_abij +einsum_permute('ijab->abij', z1);
    X2A_abij = X2A_abij +einsum_permute('iabj->abij', z2);
    X2A_abij = X2A_abij -einsum_permute('jabi->abij', z2);
    X2A_abij = X2A_abij +einsum_permute('ijab->abij', z3);
    X2A_abij = X2A_abij -einsum_permute('jiab->abij', z3);
    X2A_abij = X2A_abij +einsum_permute('bija->abij', z4);
    X2A_abij = X2A_abij -einsum_permute('aijb->abij', z4);
    X2A_abij = X2A_abij +0.5*einsum_permute('ijab->abij', z5);
    X2A_abij = X2A_abij +einsum_permute('ibja->abij', z6);
    X2A_abij = X2A_abij -einsum_permute('iajb->abij', z6);
    X2A_abij = X2A_abij -einsum_permute('jbia->abij', z6);
    X2A_abij = X2A_abij +einsum_permute('jaib->abij', z6);
    X2A_abij = X2A_abij +0.5*einsum_permute('abij->abij', z7);
    X2A_abij = X2A_abij -einsum_permute('ibja->abij', z8);
    X2A_abij = X2A_abij +einsum_permute('iajb->abij', z8);
    X2A_abij = X2A_abij +einsum_permute('jbia->abij', z8);
    X2A_abij = X2A_abij -einsum_permute('jaib->abij', z8);
    X2A_abij = X2A_abij +einsum_permute('ijab->abij', z10);
    X2A_abij = X2A_abij -einsum_permute('ijba->abij', z12);
    X2A_abij = X2A_abij +einsum_permute('ijab->abij', z12);
    X2A_abij = X2A_abij +einsum_permute('jiba->abij', z12);
    X2A_abij = X2A_abij -einsum_permute('jiab->abij', z12);
    X2A_abij = X2A_abij +einsum_permute('iabj->abij', z14);
    X2A_abij = X2A_abij +0.5*einsum_permute('ijab->abij', z18);
    X2A_abij = X2A_abij -0.5*einsum_permute('jiab->abij', z18);
    X2A_abij = X2A_abij -einsum_permute('jaib->abij', z38);
    X2A_abij = X2A_abij +einsum_permute('iajb->abij', z38);
    X2A_abij = X2A_abij -einsum_permute('jaib->abij', z44);
    X2A_abij = X2A_abij +einsum_permute('iajb->abij', z44);
    X2A_abij = X2A_abij +einsum_permute('ijab->abij', z46);
    X2A_abij = X2A_abij -einsum_permute('jiab->abij', z46);
    
    % CCSDT part
    H1A = HBar_t.H1A;
    H1B = HBar_t.H1B;
    H2A = HBar_t.H2A;
    H2B = HBar_t.H2B;
    
    D1 = einsum_kg(H1A.ov,t3a,'me,abeijm->abij');
    
    D2 = einsum_kg(H1B.ov,t3b,'me,abeijm->abij');
    
    D3 = -einsum_kg(H2B.ooov,t3b,'mnif,abfmjn->abij');
    
    D4 = -0.5*einsum_kg(H2A.ooov,t3a,'mnif,abfmjn->abij');
    
    D5 = 0.5*einsum_kg(H2A.vovv,t3a,'anef,ebfijn->abij');
    
    D6 = einsum_kg(H2B.vovv,t3b,'anef,ebfijn->abij');

    D34 = D3 + D4;
    D34 = D34 - permute(D34,[1,2,4,3]);

    D56 = D5 + D6;
    D56 = D56 - permute(D56,[2,1,3,4]);

    
    X2A_abij = X2A_abij + D1 + D2 + D34 + D56;

    for a = 1:sys.Nvir_alpha
        for b = a+1:sys.Nvir_alpha
            for i = 1:sys.Nocc_alpha
                for j = i+1:sys.Nocc_alpha
                    denom = (sys.fa_oo(i,i)+sys.fa_oo(j,j)-sys.fa_vv(a,a)-sys.fa_vv(b,b)-shift);
                    coef = X2A_abij(a,b,i,j) / denom;
                    t2a(a,b,i,j) = t2a(a,b,i,j) + coef;              
                    t2a(b,a,i,j) = -t2a(a,b,i,j);
                    t2a(a,b,j,i) = -t2a(a,b,i,j);
                    t2a(b,a,j,i) = t2a(a,b,i,j);
                end
            end
        end
    end

end