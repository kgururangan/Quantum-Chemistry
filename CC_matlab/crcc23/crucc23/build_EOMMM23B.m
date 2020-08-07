function [EOMMM23B] = build_EOMMM23B(cc_t,HBar_t,iroot)

    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C; 
    H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3D = HBar_t.H3D; H3C = HBar_t.H3C;
    VA = HBar_t.H2A.oovv; VC = HBar_t.H2C.oovv; VB = HBar_t.H2B.oovv;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;

    r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
    r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};
    
    fprintf('EOM-MMCC(2,3)B construction... ')
    
    tic
    
    D1 = einsum_kg(H3B.vvvoov,r1b,'abcije,ek->abcijk');

    D2 = einsum_kg(H3B.vvvovo,r1a,'abciek,ej->abcijk');
    D2 = D2 - permute(D2,[1,2,3,5,4,6]);

    D3 = -einsum_kg(H3B.vvoooo,r1b,'abmijk,cm->abcijk');

    D4 = -einsum_kg(H3B.vovooo,r1a,'amcijk,bm->abcijk');
    D4 = D4 - permute(D4,[2,1,3,4,5,6]);

    D5 = 0.5*einsum_kg(H3B.oovooo,r2a,'mncijk,abmn->abcijk');

    D6 = 0.5*einsum_kg(H3B.vvvvvo,r2a,'abcefk,efij->abcijk');

    D7 = einsum_kg(H3B.vooooo,r2b,'bmnjik,acmn->abcijk');
    D7 = D7 - permute(D7,[2,1,3,4,5,6]);

    D8 = einsum_kg(H3B.vvvovv,r2b,'bacjef,efik->abcijk');
    D8 = D8 - permute(D8,[1,2,3,5,4,6]);

    D9 = einsum_kg(H3A.vvooov,r2b,'abmije,ecmk->abcijk');

    D10 = einsum_kg(H3B.vvooov,r2c,'abmije,ecmk->abcijk');

    D11 = einsum_kg(H3B.vovovo,r2a,'amciek,bejm->abcijk');
    D11 = D11 - permute(D11,[2,1,3,4,5,6]) - permute(D11,[1,2,3,5,4,6]) + permute(D11,[2,1,3,5,4,6]);

    D12 = einsum_kg(H3C.vovovo,r2b,'amciek,bejm->abcijk');
    D12 = D12 - permute(D12,[2,1,3,4,5,6]) - permute(D12,[1,2,3,5,4,6]) + permute(D12,[2,1,3,5,4,6]);

    D13 = -einsum_kg(H3B.vvoovo,r2b,'banjek,ecin->abcijk');
    D13 = D13 - permute(D13,[1,2,3,5,4,6]);

    D14 = -einsum_kg(H3B.vovoov,r2b,'bmcjie,aemk->abcijk');
    D14 = D14 - permute(D14,[2,1,3,4,5,6]);

    D15 = -einsum_kg(H2B.ovoo,r2a,'mcjk,abim->abcijk');
    D15 = D15 - permute(D15,[1,2,3,5,4,6]);

    D16 = -einsum_kg(H2A.vooo,r2b,'amij,bcmk->abcijk');
    D16 = D16 - permute(D16,[2,1,3,4,5,6]);

    D17 = -einsum_kg(H2B.vooo,r2b,'amik,bcjm->abcijk');
    D17 = D17 - permute(D17,[2,1,3,4,5,6]) - permute(D17,[1,2,3,5,4,6]) + permute(D17,[2,1,3,5,4,6]);

    D18 = einsum_kg(H2B.vvvo,r2a,'bcek,aeij->abcijk');
    D18 = D18 - permute(D18,[2,1,3,4,5,6]);

    D19 = einsum_kg(H2A.vvov,r2b,'abie,ecjk->abcijk');
    D19 = D19 - permute(D19,[1,2,3,5,4,6]);

    D20 = einsum_kg(H2B.vvov,r2b,'acie,bejk->abcijk');
    D20 = D20 - permute(D20,[2,1,3,4,5,6]) - permute(D20,[1,2,3,5,4,6]) + permute(D20,[2,1,3,5,4,6]);

    I1A = einsum_kg(VA,r1a,'mnef,fn->me') + einsum_kg(VB,r1b,'mnef,fn->me');
    I1B = einsum_kg(VC,r1b,'nmfe,fn->me') + einsum_kg(VB,r1a,'nmfe,fn->me');

    D21 = -einsum_kg(einsum_kg(I1A,t2a,'me,aeij->amij'),t2b,'amij,bcmk->abcijk');
    D21 = D21 - permute(D21,[2,1,3,4,5,6]);

    D22 = -einsum_kg(einsum_kg(I1A,t2a,'me,abim->abie'),t2b,'abie,ecjk->abcijk');
    D22 = D22 - permute(D22,[1,2,3,5,4,6]);

    D23 = -einsum_kg(einsum_kg(I1B,t2b,'me,aeik->amik'),t2b,'amik,bcjm->abcijk');
    D23 = D23 - permute(D23,[2,1,3,4,5,6]) - permute(D23,[1,2,3,5,4,6]) + permute(D23,[2,1,3,5,4,6]);

    EOMMM23B = D1 + D2 + D3 + D4;

    EOMMM23B = EOMMM23B + D5 + D6 + D7 + D8;
    
    EOMMM23B = EOMMM23B + D9 + D10 + D11 + D12 + D13 + D14;

    EOMMM23B = EOMMM23B + D15 + D16 + D17 + D18 + D19 + D20;

    EOMMM23B = EOMMM23B + D21 + D22 + D23;
     
    fprintf(' finished in %4.2f s\n',toc)
     
end