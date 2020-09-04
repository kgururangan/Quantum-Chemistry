function [X2B_abij] = build_t2b_uccsdt(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, sys)

    [chi1A, chi1B, chi2A, chi2B, chi2C] = build_ucc_intermediates_v2(t1a, t1b, t2a, t2b, t2c, sys);
    [HBar_t, VT3_t] = build_ucc_hbar_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);

    % CCSD part
    X2B_abij = sys.vB_vvoo ... % 1
               -einsum_kg(chi1A.mi,t2b,'mi,abmj->abij') ... % 8
               -einsum_kg(chi1B.mj,t2b,'mj,abim->abij')... % 8
               +einsum_kg(chi1A.ae,t2b,'ae,ebij->abij')... % 8
               +einsum_kg(chi1B.be,t2b,'be,aeij->abij')... % 8
               +einsum_kg(chi2A.amie,t2b,'amie,ebmj->abij')... % 6
               +einsum_kg(chi2B.amie,t2c,'amie,ebmj->abij')... % 6
               +einsum_kg(chi2B.mbej,t2a,'mbej,aeim->abij')... % 4
               +einsum_kg(chi2C.bmje,t2b,'bmje,aeim->abij')... % 4
               +einsum_kg(chi2B.mnij,t2b,'mnij,abmn->abij')... % 5
               +einsum_kg(chi2B.abef,t2b,'abef,efij->abij')... % 4
               -einsum_kg(chi2B.amej,t2b,'amej,ebim->abij')... % 4
               -einsum_kg(chi2B.mbie,t2b,'mbie,aemj->abij')... % 5
               +einsum_kg(chi2B.abej,t1a,'abej,ei->abij')... % 7
               +einsum_kg(chi2B.abie,t1b,'abie,ej->abij')... % 5
               -einsum_kg(chi2B.mbij,t1a,'mbij,am->abij')... % 2
               -einsum_kg(sys.vB_vooo,t1b,'amij,bm->abij'); % 1


    % CCSDT part
    H1A = HBar_t.H1A;
    H1B = HBar_t.H1B;
    H2A = HBar_t.H2A;
    H2B = HBar_t.H2B;
    H2C = HBar_t.H2C;
    
    D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
    D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
    D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
    D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
    D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
    D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
    D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
    D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
    D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
    D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
      
    X2B_abij = X2B_abij + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10;

end