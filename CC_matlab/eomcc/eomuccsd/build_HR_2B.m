function [X2B] = build_HR_2B(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys)


    Nunocc_a = sys.Nvir_alpha; Nocc_a = sys.Nocc_alpha;
    Nunocc_b = sys.Nvir_beta; Nocc_b = sys.Nocc_beta;

    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

    t2b = cc_t.t2b; 

    X2B = zeros(Nunocc_a, Nunocc_b, Nocc_a, Nocc_b);

    X2B = X2B + einsum_kg(H1A.vv,r2b,'ae,ebij->abij');
    X2B = X2B + einsum_kg(H1B.vv,r2b,'be,aeij->abij');
    X2B = X2B - einsum_kg(H1A.oo,r2b,'mi,abmj->abij');
    X2B = X2B - einsum_kg(H1B.oo,r2b,'mj,abim->abij');

    X2B = X2B + einsum_kg(H2B.oooo,r2b,'mnij,abmn->abij');
    X2B = X2B + einsum_kg(H2B.vvvv,r2b,'abef,efij->abij');

    X2B = X2B + einsum_kg(H2A.voov,r2b,'amie,ebmj->abij');
    X2B = X2B + einsum_kg(H2B.voov,r2c,'amie,ebmj->abij');
    X2B = X2B + einsum_kg(H2B.ovvo,r2a,'mbej,aeim->abij');
    X2B = X2B + einsum_kg(H2C.voov,r2b,'bmje,aeim->abij');
    X2B = X2B - einsum_kg(H2B.ovov,r2b,'mbie,aemj->abij');
    X2B = X2B - einsum_kg(H2B.vovo,r2b,'amej,ebim->abij');
% 
    X2B = X2B + einsum_kg(H2B.vvvo,r1a,'abej,ei->abij');
    X2B = X2B + einsum_kg(H2B.vvov,r1b,'abie,ej->abij');
    X2B = X2B - einsum_kg(H2B.ovoo,r1a,'mbij,am->abij');
    X2B = X2B - einsum_kg(H2B.vooo,r1b,'amij,bm->abij');

%     X2B = X2B - 0.5*einsum_kg(H3B.oovovo,r2a,'mnbifj,afmn->abij');
%     X2B = X2B - einsum_kg(H3C.ovooov,r2b,'mbnijf,afmn->abij');
%     X2B = X2B - einsum_kg(H3B.vooovo,r2b,'anmifj,fbnm->abij');
%     X2B = X2B - 0.5*einsum_kg(H3C.voooov,r2c,'amnijf,fbnm->abij');
% 
%     X2B = X2B + 0.5*einsum_kg(H3B.vovvvo,r2a,'anbefj,efin->abij');
%     X2B = X2B + einsum_kg(H3C.vvovov,r2b,'abnejf,efin->abij');
%     X2B = X2B + einsum_kg(H3B.vovovv,r2b,'anbife,fenj->abij');
%     X2B = X2B + 0.5*einsum_kg(H3C.vvoovv,r2c,'abnief,efjn->abij');

    X2B = X2B - 0.5*einsum_kg(einsum_kg(sys.vA_oovv,r2a,'mnef,afmn->ae'),t2b,'ae,ebij->abij')...
              - 0.5*einsum_kg(einsum_kg(sys.vA_oovv,r2a,'mnef,efin->mi'),t2b,'mi,abmj->abij');

    X2B = X2B - einsum_kg(einsum_kg(sys.vB_oovv,r2b,'nmfe,fbnm->be'),t2b,'be,aeij->abij')...
              - einsum_kg(einsum_kg(sys.vB_oovv,r2b,'mnef,afmn->ae'),t2b,'ae,ebij->abij')...
              - einsum_kg(einsum_kg(sys.vB_oovv,r2b,'nmfe,fenj->mj'),t2b,'mj,abim->abij')...
              - einsum_kg(einsum_kg(sys.vB_oovv,r2b,'mnef,efin->mi'),t2b,'mi,abmj->abij');

    X2B = X2B - 0.5*einsum_kg(einsum_kg(sys.vC_oovv,r2c,'mnef,bfmn->be'),t2b,'be,aeij->abij')...
              - 0.5*einsum_kg(einsum_kg(sys.vC_oovv,r2c,'mnef,efjn->mj'),t2b,'mj,abim->abij');

%     X2B = X2B + einsum_kg(H3B.vovovo,r1a,'ambiej,em->abij');
%     X2B = X2B + einsum_kg(H3C.vvooov,r1b,'abmije,em->abij');

    X2B = X2B + einsum_kg(einsum_kg(H2B.ovvv,r1a,'mbef,em->bf'),t2b,'bf,afij->abij')...
              - einsum_kg(einsum_kg(H2B.oovo,r1a,'mnej,em->nj'),t2b,'nj,abin->abij')...
              + einsum_kg(einsum_kg(H2A.vovv,r1a,'amfe,em->af'),t2b,'af,fbij->abij')...
              - einsum_kg(einsum_kg(H2A.ooov,r1a,'nmie,em->ni'),t2b,'ni,abnj->abij');

    X2B = X2B + einsum_kg(einsum_kg(H2B.vovv,r1b,'amfe,em->af'),t2b,'af,fbij->abij')...
              - einsum_kg(einsum_kg(H2B.ooov,r1b,'nmie,em->ni'),t2b,'ni,abnj->abij')...
              + einsum_kg(einsum_kg(H2C.vovv,r1b,'bmfe,em->bf'),t2b,'bf,afij->abij')...
              - einsum_kg(einsum_kg(H2C.oovo,r1b,'mnej,em->nj'),t2b,'nj,abin->abij');

end

