function [t2b] = update_t2b_jun(t1a,t1b,t2a,t2b,t2c,sys)

    s17 = einsum_kg(sys.vB_ovvo,t1a,'mbej,ei->jbmi');
    s19 = einsum_kg(sys.vB_oooo,t1a,'mnij,am->ijna');
    s21 = einsum_kg(sys.vB_ovov,t1b,'mbie,ej->ibmj');
    s23 = einsum_kg(sys.vB_vovo,t1a,'amej,ei->jami');
    s25 = einsum_kg(sys.vB_vvvv,t1a,'abef,ei->abfi');
    q1 = einsum_kg(sys.fa_ov,t1a,'me,ei->mi');
    q2 = einsum_kg(sys.fa_ov,t1a,'me,am->ea');
    s29 = einsum_kg(sys.vA_ooov,t2b,'mnie,ebnj->imjb');
    q3 = einsum_kg(sys.vA_ooov,t1a,'mnie,en->im');
    s32 = einsum_kg(sys.vA_ovvv,t1a,'maef,ei->afmi');
    q4 = einsum_kg(sys.vA_ovvv,t1a,'maef,fm->ae');
    s35 = einsum_kg(sys.vB_oovo,t1a,'mnej,ei->jmni');
    s93 = einsum_kg(s35,t1a,'jmni,am->ijna');
    s39 = einsum_kg(sys.vB_oovo,t2b,'mnej,ebin->jmib');
    q5 = einsum_kg(sys.vB_oovo,t1a,'mnej,em->jn');
    s42 = einsum_kg(sys.vB_ovvv,t2b,'mbef,efij->bmij');
    q6 = einsum_kg(sys.vB_ovvv,t1a,'mbef,em->bf');
    s45 = einsum_kg(sys.vB_ooov,t2c,'mnie,bejn->imjb');
    s47 = einsum_kg(sys.vB_vovv,t1a,'amef,ei->afmi');
    s107 = einsum_kg(s47,t1b,'afmi,fj->iamj');
    s51 = einsum_kg(sys.vB_oovo,t2a,'mnej,aeim->jnia');
    s53 = einsum_kg(sys.vB_ovvv,t1b,'mbef,fj->bemj');
    q7 = einsum_kg(sys.fb_ov,t1b,'me,ej->mj');
    q8 = einsum_kg(sys.fb_ov,t1b,'me,bm->eb');
    s57 = einsum_kg(sys.vB_ooov,t1b,'mnie,ej->imnj');
    s105 = einsum_kg(s57,t1a,'imnj,am->ijna');
    q9 = einsum_kg(sys.vB_ooov,t1b,'mnie,en->im');
    s62 = einsum_kg(sys.vB_vovv,t1b,'amef,fj->aemj');
    s64 = einsum_kg(sys.vB_vovv,t2b,'amef,efij->amij');
    q10 = einsum_kg(sys.vB_vovv,t1b,'amef,fm->ae');
    s67 = einsum_kg(sys.vA_ooov,t2b,'mnje,aeim->jnia');
    q11 = einsum_kg(sys.vA_ooov,t1b,'mnje,en->jm');
    s70 = einsum_kg(sys.vA_ovvv,t1b,'mbef,fj->bemj');
    q12 = einsum_kg(sys.vA_ovvv,t1b,'mbef,fm->be');
    s73 = einsum_kg(sys.vA_oovv,t2a,'mnef,aeim->fnia');
    q13 = einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi');
    q14 = einsum_kg(sys.vA_oovv,t2a,'mnef,aemn->fa');
    s77 = einsum_kg(sys.vB_oovv,t2a,'mnef,aeim->fnia');
    s122 = einsum_kg(s77,t1b,'fnia,fj->ianj');
    q16 = einsum_kg(sys.vB_oovv,t2b,'mnef,ebmn->fb');
    s81 = einsum_kg(sys.vB_oovv,t2b,'mnef,efij->mnij');
    s118 = einsum_kg(s81,t1a,'mnij,am->ijna');
    s84 = einsum_kg(sys.vB_oovv,t2b,'mnef,afmj->enja');
    s111 = einsum_kg(s84,t1a,'enja,ei->jani');
    q18 = einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi');
    q19 = einsum_kg(sys.vA_oovv,t2c,'mnef,bfmn->eb');
    q20 = einsum_kg(sys.vA_oovv,t2c,'mnef,efjn->mj');
    s91 = einsum_kg(sys.vA_oovv,t2b,'mnef,aeim->fnia');
    s124 = einsum_kg(s91,t1b,'fnia,fj->ianj');
    s37 = einsum_kg(sys.vB_ovvv,t1a,'mbef,ei->bfmi');
    s95 = einsum_kg(s37,t1b,'bfmi,fj->ibmj');
    s97 = einsum_kg(sys.vA_oovv,t1a,'mnef,ei->fmni');
    s98 = einsum_kg(s97,t2b,'fmni,fbnj->imjb');
    q21 = einsum_kg(s97,t1a,'fmni,fn->im');
    s102 = einsum_kg(sys.vB_oovv,t1a,'mnef,ei->fmni');
    s109 = einsum_kg(s102,t1b,'fmni,fj->imnj');
    q24 = einsum_kg(s102,t1b,'fmni,fn->im');
    s103 = einsum_kg(s102,t2c,'fmni,bfjn->imjb');
    s128 = einsum_kg(s109,t1a,'imnj,am->ijna');
    s86 = einsum_kg(sys.vB_oovv,t2b,'mnef,afin->emia');
    s49 = einsum_kg(sys.vB_ovvo,t1b,'ieam,ej->iamj');
    s114 = einsum_kg(sys.vB_oovv,t1b,'mnef,fj->emnj');
    q25 = einsum_kg(s114,t1a,'emnj,em->jn');
    s115 = einsum_kg(s114,t2b,'emnj,ebin->jmib');
    q17 = einsum_kg(sys.vB_oovv,t2b,'mnef,efmj->nj');
    q26 = einsum_kg(sys.vB_oovv,t1a,'mnef,em->fn');
    q27 = einsum_kg(q26,t1b,'fn,bn->fb');
    q28 = einsum_kg(sys.vB_oovv,t1b,'mnef,fn->em');
    q29 = einsum_kg(q28,t1a,'em,am->ea');
    q15 = einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ea');
    q22 = einsum_kg(sys.vA_oovv,t1a,'mnef,fn->em');
    q23 = einsum_kg(q22,t1a,'em,am->ea');
    q30 = einsum_kg(sys.vA_oovv,t1b,'mnef,fn->em');
    q32 = einsum_kg(q30,t1b,'em,bm->eb');
    q31 = einsum_kg(q30,t1b,'em,ej->mj');
    s59 = einsum_kg(sys.vB_ooov,t2b,'mnie,aemj->inja');


    x1 = einsum_permute(sys.vB_ooov,'ijmb->ijbm') + ...
         einsum_permute('jbmi->ijbm', s17) + ...
         einsum_permute('ibmj->ijbm', s21) + ...
         einsum_permute('imjb->ijbm', s29) - ...
         einsum_permute('jmib->ijbm', s39) + ...
         einsum_permute('bmij->ijbm', s42) + ...
         einsum_permute('imjb->ijbm', s45) + ...
         einsum_permute('ibmj->ijbm', s95) + ...
         einsum_permute('imjb->ijbm', s98) + ...
         einsum_permute('imjb->ijbm', s103) - ...
         einsum_permute('jmib->ijbm', s115);

    z1 = einsum_kg(x1,t1a,'ijbm,am->ijba');
    z2 = einsum_kg(sys.vB_vovv,t1a,'ejab,ei->jabi');

    x2 = einsum_permute('ijam->ijam', sys.vB_oovo) - ...
         einsum_permute('ijma->ijam', s19) + ...
         einsum_permute('jami->ijam', s23) - ...
         einsum_permute('ijma->ijam', s93) + ...
         einsum_permute('iamj->ijam', s107) + ...
         einsum_permute('jmia->ijam', s51) - ...
         einsum_permute('ijma->ijam', s105) + ...
         einsum_permute('amij->ijam', s64) - ...
         einsum_permute('jmia->ijam', s67) + ...
         einsum_permute('iamj->ijam', s122) - ...
         einsum_permute('ijma->ijam', s118) - ...
         einsum_permute('jami->ijam', s111) + ...
         einsum_permute('iamj->ijam', s124) - ...
         einsum_permute('ijma->ijam', s128) + ...
         einsum_permute('iamj->ijam', s49) - ...
         einsum_permute('imja->ijam', s59);

    z3 = einsum_kg(x2,t1b,'ijam,bm->ijab');

    x3 = einsum_permute(sys.vB_ovvv,'ieab->iabe') + einsum_permute('abei->iabe', s25);

    z4 = einsum_kg(x3,t1b,'iabe,ej->iabj');

    x4 = einsum_permute(sys.vB_ovvo,'mbej->jbem') + einsum_permute('bemj->jbem', s53);

    z5 = einsum_kg(x4,t2a,'jbem,aeim->jbia');

    x5 = einsum_permute(sys.fa_oo,'mi->im') + ...
         einsum_permute('mi->im', q1) + ...
         einsum_permute('im->im', q3) + ...
         einsum_permute('im->im', q9) + ...
         0.5*einsum_permute('mi->im', q13) + ...
         einsum_permute('mi->im', q18) + ...
         einsum_permute('im->im', q21) + ...
         einsum_permute('im->im', q24);

    z6 = einsum_kg(x5,t2b,'im,abmj->ijab');

    x6 = einsum_permute(sys.fa_vv,'ae->ae') - ...
         einsum_permute('ea->ae', q2) - ...
         einsum_permute('ae->ae', q4) + ...
         einsum_permute('ae->ae', q10) + ...
         0.5*einsum_permute('ea->ae', q14) - ...
         einsum_permute('ea->ae', q29) - ...
         einsum_permute('ea->ae', q15) - ...
         einsum_permute('ea->ae', q23);

    z7 = einsum_kg(x6,t2b,'ae,ebij->aijb');

    x7 = einsum_permute(sys.fb_oo,'mj->jm') + ...
         einsum_permute('jm->jm', q5) + ...
         einsum_permute('mj->jm', q7) + ...
         einsum_permute('jm->jm', q11) + ...
         0.5*einsum_permute('mj->jm', q20) + ...
         einsum_permute('jm->jm', q25) + ...
         einsum_permute('mj->jm', q17) + ...
         einsum_permute('mj->jm', q31);

    z8 = einsum_kg(x7,t2b,'jm,abim->jiab');

    x8 = einsum_permute(sys.fb_vv,'be->be') + ...
         einsum_permute('be->be', q6) - ...
         einsum_permute('eb->be', q8) - ...
         einsum_permute('be->be', q12) - ...
         einsum_permute('eb->be', q16) - ...
         0.5*einsum_permute('eb->be', q19) - ...
         einsum_permute('eb->be', q27) - ...
         einsum_permute('eb->be', q32);

    z9 = einsum_kg(x8,t2b,'be,aeij->bija');

    x9 = einsum_permute(sys.vA_ovov,'maie->iaem') + ...
         einsum_permute('aemi->iaem', s32) - ...
         einsum_permute('emia->iaem', s73) - ...
         einsum_permute('emia->iaem', s86);

    z10 = einsum_kg(x9,t2b,'iaem,ebmj->iajb');

    x10 = einsum_permute(sys.vB_oooo,'mnij->ijmn') + ...
          einsum_permute('jmni->ijmn', s35) + ...
          einsum_permute('imnj->ijmn', s57) + ...
          einsum_permute('mnij->ijmn', s81) + ...
          einsum_permute('imnj->ijmn', s109);

    z11 = einsum_kg(x10,t2b,'ijmn,abmn->ijab');
    x11 = einsum_permute('mbie->ibem', sys.vB_ovov) + einsum_permute('bemi->ibem', s37);

    z12 = einsum_kg(x11,t2b,'ibem,aemj->ibja');

    x12 = einsum_permute(sys.vB_vovo,'amej->jaem') + einsum_permute('aemj->jaem', s62) - einsum_permute('emja->jaem', s84);

    z13 = einsum_kg(x12,t2b,'jaem,ebim->jaib');

    z14 = einsum_kg(sys.vB_vvvv,t2b,'abef,efij->abij');

    x13 = einsum_permute(sys.vA_ovov,'mbje->jbem') - einsum_permute('bemj->jbem', s70);

    z15 = einsum_kg(x13,t2b,'jbem,aeim->jbia');

    x14 = einsum_permute(sys.vB_ovvo,'ieam->iaem') + einsum_permute('aemi->iaem', s47) + einsum_permute('emia->iaem', s77) + einsum_permute('emia->iaem', s91);

    z16 = einsum_kg(x14,t2c,'iaem,bejm->iajb');


    new_t = einsum_permute(sys.vB_vvoo,'abij->abij');
    new_t = new_t - einsum_permute(z1,'ijba->abij');
    new_t = new_t + einsum_permute(z2,'jabi->abij');
    new_t = new_t - einsum_permute(z3,'ijab->abij');
    new_t = new_t + einsum_permute(z4,'iabj->abij');
    new_t = new_t + einsum_permute(z5,'jbia->abij');
    new_t = new_t - einsum_permute(z6,'ijab->abij');
    new_t = new_t + einsum_permute('aijb->abij', z7);
    new_t = new_t - einsum_permute('jiab->abij', z8);
    new_t = new_t + einsum_permute('bija->abij', z9);
    new_t = new_t - einsum_permute('iajb->abij', z10);
    new_t = new_t + einsum_permute('ijab->abij', z11);
    new_t = new_t - einsum_permute('ibja->abij', z12);
    new_t = new_t - einsum_permute('jaib->abij', z13);
    new_t = new_t + einsum_permute('abij->abij', z14);
    new_t = new_t - einsum_permute('jbia->abij', z15);
    new_t = new_t + einsum_permute('iajb->abij', z16);

    for i = 1:sys.Nocc_alpha
        for j = 1:sys.Nocc_beta
            for a = 1:sys.Nvir_alpha
                for b = 1:sys.Nvir_beta
                    denom = (sys.fa_oo(i,i)+sys.fb_oo(j,j)-sys.fa_vv(a,a)-sys.fb_vv(b,b)); 
                    %coef = -new_t(a,b,i,j) / denom;
                    coef = new_t(a,b,i,j) / denom;
                    t2b(a,b,i,j) = t2b(a,b,i,j) + coef;
                end
            end
        end
    end



end
