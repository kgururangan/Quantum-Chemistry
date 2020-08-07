function [EOMMM23C] = build_EOMMM23C(cc_t,HBar_t,iroot)

    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C; 
    H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3D = HBar_t.H3D; H3C = HBar_t.H3C;
    VA = HBar_t.H2A.oovv; VC = HBar_t.H2C.oovv; VB = HBar_t.H2B.oovv;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;

    r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
    r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};
    
    fprintf('EOM-MMCC(2,3)C construction... ')
    
    tic

    
    D1 = -einsum_kg(H3C.ovvooo,r1a,'mbcijk,am->abcijk');

    D2 = -einsum_kg(H3C.vvoooo,r1b,'acmikj,bm->abcijk');
    D2 = D2 - permute(D2,[1,3,2,4,5,6]);

    D3 = einsum_kg(H3C.vvvvoo,r1a,'abcejk,ei->abcijk');
    
    D4 = einsum_kg(H3C.vvvoov,r1b,'acbike,ej->abcijk');
    D4 = D4 - permute(D4,[1,2,3,4,6,5]);

    D5 = 0.5*einsum_kg(H3C.vooooo,r2c,'amnijk,bcmn->abcijk');

    D6 = einsum_kg(H3C.oovooo,r2b,'mncijk,abmn->abcijk');
    D6 = D6 - permute(D6,[1,3,2,4,5,6]);

    D7 = 0.5*einsum_kg(H3C.vvvovv,r2c,'abcief,efjk->abcijk');
    
    D8 = einsum_kg(H3C.vvvvvo,r2b,'abcefk,efij->abcijk');
    D8 = D8 - permute(D8,[1,2,3,4,6,5]);


    D9 = einsum_kg(H3C.ovvvoo,r2a,'mbcejk,aeim->abcijk');

    D10 = einsum_kg(H3D.vovovo,r2b,'bmcjek,aeim->abcijk');
    
    D11 = einsum_kg(H3B.vovovo,r2b,'ambiej,ecmk->abcijk');
    D11 = D11 - permute(D11,[1,3,2,4,5,6]) - permute(D11,[1,2,3,4,6,5]) + permute(D11,[1,3,2,4,6,5]);

    D12 = einsum_kg(H3C.vovovo,r2c,'ambiej,ecmk->abcijk');
    D12 = D12 - permute(D12,[1,3,2,4,5,6]) - permute(D12,[1,2,3,4,6,5]) + permute(D12,[1,3,2,4,6,5]);

    D13 = -einsum_kg(H3C.ovvovo,r2b,'mbciek,aemj->abcijk');
    D13 = D13 - permute(D13,[1,2,3,4,6,5]);

    D14 = -einsum_kg(H3C.vvovoo,r2b,'acmekj,ebim->abcijk');
    D14 = D14 - permute(D14,[1,3,2,4,5,6]);

    D15 = -einsum_kg(H2B.ovoo,r2b,'mcik,abmj->abcijk');
    D15 = D15 - permute(D15,[1,3,2,4,5,6]) - permute(D15,[1,2,3,4,6,5]) + permute(D15,[1,3,2,4,6,5]);

    D16 = -einsum_kg(H2B.vooo,r2c,'amij,bcmk->abcijk');
    D16 = D16 - permute(D16,[1,2,3,4,6,5]);

    D17 = -einsum_kg(H2C.vooo,r2b,'cmkj,abim->abcijk');
    D17 = D17 - permute(D17,[1,3,2,4,5,6]);

    D18 = einsum_kg(H2B.vvvo,r2b,'acek,ebij->abcijk');
    D18 = D18 - permute(D18,[1,3,2,4,5,6]) - permute(D18,[1,2,3,4,6,5]) + permute(D18,[1,3,2,4,6,5]);

    D19 = einsum_kg(H2B.vvov,r2c,'abie,ecjk->abcijk');
    D19 = D19 - permute(D19,[1,3,2,4,5,6]);

    D20 = einsum_kg(H2C.vvov,r2b,'cbke,aeij->abcijk');
    D20 = D20 - permute(D20,[1,2,3,4,6,5]);


    I1A = einsum_kg(VA,r1a,'mnef,fn->me') + einsum_kg(VB,r1b,'mnef,fn->me');
    I1B = einsum_kg(VC,r1b,'nmfe,fn->me') + einsum_kg(VB,r1a,'nmfe,fn->me');

    D21 = -einsum_kg(einsum_kg(I1B,t2c,'me,ecjk->mcjk'),t2b,'mcjk,abim->abcijk');
    D21 = D21 - permute(D21,[1,3,2,4,5,6]);

    D22 = -einsum_kg(einsum_kg(I1B,t2c,'me,bcmk->bcek'),t2b,'bcek,aeij->abcijk');
    D22 = D22 - permute(D22,[1,2,3,4,6,5]);

    D23 = -einsum_kg(einsum_kg(I1A,t2b,'me,ecik->mcik'),t2b,'mcik,abmj->abcijk');
    D23 = D23 - permute(D23,[1,3,2,4,5,6]) - permute(D23,[1,2,3,4,6,5]) + permute(D23,[1,3,2,4,6,5]);

     EOMMM23C = D1 + D2 + D3 + D4;
% 
     EOMMM23C = EOMMM23C + D5 + D6 + D7 + D8;
%     
     EOMMM23C = EOMMM23C + D9 + D10 + D11 + D12 + D13 + D14;
% 
     EOMMM23C = EOMMM23C + D15 + D16 + D17 + D18 + D19 + D20;

    EOMMM23C = EOMMM23C + D21 + D22 + D23;
     
    fprintf(' finished in %4.2f s\n',toc)
     
end