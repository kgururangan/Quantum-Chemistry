function [EOMMM23C] = build_EOMMM23C(cc_t,HBar_t,iroot)

    %H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C; 
    %H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3D = HBar_t.H3D; H3C = HBar_t.H3C;
    VA = HBar_t.H2A.oovv; VC = HBar_t.H2C.oovv; VB = HBar_t.H2B.oovv;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;

    r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
    r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};
    
    fprintf('EOM-MMCC(2,3)C construction... ')
    
    tic
   
    %D1 = -einsum_kg(H3C.ovvooo,r1a,'mbcijk,am->abcijk');
    I1 = -einsum_kg(H2B.ovov,r1a,'mbie,am->abie');
    I2 = -einsum_kg(H2B.oooo,r1a,'mnij,am->anij');
    I3 = -einsum_kg(H2B.ovvo,r1a,'mbej,am->abej');
    A1 = einsum_kg(I1,t2c,'abie,ecjk->abcijk');
    A2 = -einsum_kg(I2,t2c,'anij,bcnk->abcijk');
    A3 = einsum_kg(I3,t2b,'abej,ecik->abcijk');
    A1 = A1 - permute(A1,[1,3,2,4,5,6]);
    A2 = A2 - permute(A2,[1,2,3,4,6,5]);
    A3 = A3 - permute(A3,[1,3,2,4,5,6]) - permute(A3,[1,2,3,4,6,5]) + permute(A3,[1,3,2,4,6,5]);
    D1 = A1 + A2 + A3;

    %D2 = -einsum_kg(H3C.vvoooo,r1b,'acmikj,bm->abcijk');
    %D2 = D2 - permute(D2,[1,3,2,4,5,6]);
    I1 = -einsum_kg(H2B.vovo,r1b,'amej,bm->abej');
    I2 = -einsum_kg(H2B.voov,r1b,'amie,bm->abie');
    I3 = -einsum_kg(H2C.voov,r1b,'cmke,bm->bcek'); I3 = I3 - permute(I3,[2,1,3,4]);
    I4 = -einsum_kg(H2C.oooo,r1b,'mnjk,bm->bnjk');
    I5 = -einsum_kg(H2B.oooo,r1b,'nmij,bm->nbij');
    A1 = einsum_kg(I1,t2b,'abej,ecik->abcijk');
    A2 = einsum_kg(I2,t2c,'abie,ecjk->abcijk');
    A3 = einsum_kg(I3,t2b,'bcek,aeij->abcijk');
    A4 = -einsum_kg(I4,t2b,'bnjk,acin->abcijk');
    A5 = -einsum_kg(I5,t2b,'nbij,acnk->abcijk');
    A15 = A1 + A5;
    A24 = A2 + A4;
    A3 = A3 - permute(A3,[1,2,3,4,6,5]);
    A24 = A24 - permute(A24,[1,3,2,4,5,6]);
    A15 = A15 - permute(A15,[1,3,2,4,5,6]) - permute(A15,[1,2,3,4,6,5]) + permute(A15,[1,3,2,4,6,5]);
    D2 = A3 + A24 + A15;

    %D3 = einsum_kg(H3C.vvvvoo,r1a,'abcejk,ei->abcijk');
    I1 = einsum_kg(H2B.vovo,r1a,'amej,ei->amij');
    A1 = -einsum_kg(I1,t2c,'amij,bcmk->abcijk');
    A1 = A1 - permute(A1,[1,2,3,4,6,5]);
    I2 = einsum_kg(H2B.vvvv,r1a,'abef,ei->abif');
    A2 = einsum_kg(I2,t2c,'abif,fcjk->abcijk');
    A2 = A2 - permute(A2,[1,3,2,4,5,6]);
    I3 = einsum_kg(H2B.ovvo,r1a,'mbej,ei->mbij');
    A3 = -einsum_kg(I3,t2b,'mbij,acmk->abcijk');
    A3 = A3 - permute(A3,[1,3,2,4,5,6]) - permute(A3,[1,2,3,4,6,5]) + permute(A3,[1,3,2,4,6,5]);
    D3 = A1 + A2 + A3;
    
    %D4 = einsum_kg(H3C.vvvoov,r1b,'acbike,ej->abcijk');
    %D4 = D4 - permute(D4,[1,2,3,4,6,5]);
    I1 = einsum_kg(H2B.ovov,r1b,'mbie,ej->mbij');
    A1 = -einsum_kg(I1,t2b,'mbij,acmk->abcijk');
    I2 = -einsum_kg(H2C.voov,r1b,'bmke,ej->bmkj'); I2 = I2 - permute(I2,[1,2,4,3]);
    A2 = -einsum_kg(I2,t2b,'bmkj,acim->abcijk');
    I3 = einsum_kg(H2B.voov,r1b,'amie,ej->amij');
    A3 = -einsum_kg(I3,t2c,'amij,bcmk->abcijk');
    I4 = einsum_kg(H2C.vvvv,r1b,'bcef,ej->bcjf');
    A4 = einsum_kg(I4,t2b,'bcjf,afik->abcijk');
    I5 = einsum_kg(H2B.vvvv,r1b,'abfe,ej->abfj');
    A5 = einsum_kg(I5,t2b,'abfj,fcik->abcijk');
    A34 = A3 + A4;
    A15 = A1 + A5;
    A34 = A34 - permute(A34,[1,2,3,4,6,5]);
    A2 = A2 - permute(A2,[1,3,2,4,5,6]);
    A15 = A15 - permute(A15,[1,3,2,4,5,6]) - permute(A15,[1,2,3,4,6,5]) + permute(A15,[1,3,2,4,6,5]);
    D4 = A2 + A15 + A34;
    
    %D5 = 0.5*einsum_kg(H3C.vooooo,r2c,'amnijk,bcmn->abcijk');
    I1 = 0.5*einsum_kg(H2C.ooov,r2c,'mnje,bcmn->bcje');
    A1 = einsum_kg(I1,t2b,'bcje,aeik->abcijk');
    A1 = A1 - permute(A1,[1,2,3,4,6,5]);
    D5 = A1;

    %D6 = einsum_kg(H3C.oovooo,r2b,'mncijk,abmn->abcijk');
    %D6 = D6 - permute(D6,[1,3,2,4,5,6]);
    I1 = einsum_kg(H2B.ooov,r2b,'mnie,abmn->abie');
    A1 = einsum_kg(I1,t2c,'abie,ecjk->abcijk');
    A1 = A1 - permute(A1,[1,3,2,4,5,6]);
    I2 = einsum_kg(H2B.oovo,r2b,'mnej,abmn->abej');
    A2 = einsum_kg(I2,t2b,'abej,ecik->abcijk');
    A2 = A2 - permute(A2,[1,3,2,4,5,6]) - permute(A2,[1,2,3,4,6,5]) + permute(A2,[1,3,2,4,6,5]);
    D6 = A1 + A2;

    %D7 = 0.5*einsum_kg(H3C.vvvovv,r2c,'abcief,efjk->abcijk');
    I1 = 0.5*einsum_kg(H2C.ovvv,r2c,'mbfe,efjk->mbkj');
    A1 = -einsum_kg(I1,t2b,'mbkj,acim->abcijk');
    A1 = A1 - permute(A1,[1,3,2,4,5,6]);
    D7 = A1;
    
    %D8 = einsum_kg(H3C.vvvvvo,r2b,'abcefk,efij->abcijk');
    %D8 = D8 - permute(D8,[1,2,3,4,6,5]);
    I1 = einsum_kg(H2B.vovv,r2b,'amef,efij->amij');
    A1 = -einsum_kg(I1,t2c,'amij,bcmk->abcijk');
    A1 = A1 - permute(A1,[1,2,3,4,6,5]);
    I2 = einsum_kg(H2B.ovvv,r2b,'mbef,efij->mbij');
    A2 = -einsum_kg(I2,t2b,'mbij,acmk->abcijk');
    A2 = A2 - permute(A2,[1,3,2,4,5,6]) - permute(A2,[1,2,3,4,6,5]) + permute(A2,[1,3,2,4,6,5]);
    D8 = A1 + A2;

    %D9 = einsum_kg(H3C.ovvvoo,r2a,'mbcejk,aeim->abcijk');
    I1 = einsum_kg(H2B.oovo,r2a,'mnej,aeim->anij');
    A1 = -einsum_kg(I1,t2c,'anij,bcnk->abcijk');
    A1 = A1 - permute(A1,[1,2,3,4,6,5]);
    I2 = einsum_kg(H2B.ovvv,r2a,'mbef,aeim->abif');
    A2 = einsum_kg(I2,t2c,'abif,fcjk->abcijk');
    A2 = A2 - permute(A2,[1,3,2,4,5,6]);
    D9 = A1 + A2;

    %D10 = einsum_kg(H3D.vovovo,r2b,'bmcjek,aeim->abcijk');
    I1 = einsum_kg(H2C.oovo,r2b,'mnej,aeim->anij');
    A1 = -einsum_kg(I1,t2c,'anij,bcnk->abcijk');
    A1 = A1 - permute(A1,[1,2,3,4,6,5]);
    I2 = einsum_kg(H2C.vovv,r2b,'bmfe,aeim->abif');
    A2 = einsum_kg(I2,t2c,'abif,fcjk->abcijk');
    A2 = A2 - permute(A2,[1,3,2,4,5,6]);
    D10 = A1 + A2;
    
    %D11 = einsum_kg(H3B.vovovo,r2b,'ambiej,ecmk->abcijk');
    %D11 = D11 - permute(D11,[1,3,2,4,5,6]) - permute(D11,[1,2,3,4,6,5]) + permute(D11,[1,3,2,4,6,5]);
    I1 = einsum_kg(H2A.vovv,r2b,'amfe,ecmk->acfk');
    A1 = einsum_kg(I1,t2b,'acfk,fbij->abcijk');
    I2 = einsum_kg(H2B.ovvv,r2b,'mbef,ecmk->bcfk'); I2 = I2 - permute(I2,[2,1,3,4]);
    A2 = einsum_kg(I2,t2b,'bcfk,afij->abcijk');
    A2 = A2 - permute(A2,[1,2,3,4,6,5]);
    I3 = einsum_kg(H2A.oovo,r2b,'mnei,ecmk->ncik'); 
    A3 = -einsum_kg(I3,t2b,'ncik,abnj->abcijk');
    I4 = einsum_kg(H2B.oovo,r2b,'mnej,ecmk->cnkj'); I4 = I4 - permute(I4,[1,2,4,3]);
    A4 = -einsum_kg(I4,t2b,'cnkj,abin->abcijk');
    A4 = A4 - permute(A4,[1,3,2,4,5,6]);
    A13 = A1 + A3;
    A13 = A13 - permute(A13,[1,3,2,4,5,6]) - permute(A13,[1,2,3,4,6,5]) + permute(A13,[1,3,2,4,6,5]);
    D11 = A13 + A2 + A4;

    %D12 = einsum_kg(H3C.vovovo,r2c,'ambiej,ecmk->abcijk');
    %D12 = D12 - permute(D12,[1,3,2,4,5,6]) - permute(D12,[1,2,3,4,6,5]) + permute(D12,[1,3,2,4,6,5]);
    I1 = einsum_kg(H2B.vovv,r2c,'amfe,cekm->acfk');
    A1 = einsum_kg(I1,t2b,'acfk,fbij->abcijk');
    I2 = einsum_kg(H2C.vovv,r2c,'bmfe,cekm->bcfk'); I2 = I2 - permute(I2,[2,1,3,4]);
    A2 = einsum_kg(I2,t2b,'bcfk,afij->abcijk');
    A2 = A2 - permute(A2,[1,2,3,4,6,5]);
    I3 = einsum_kg(H2B.ooov,r2c,'nmie,cekm->ncik');
    A3 = -einsum_kg(I3,t2b,'ncik,abnj->abcijk');
    I4 = einsum_kg(H2C.ooov,r2c,'nmje,cekm->ncjk'); I4 = I4 - permute(I4,[1,2,4,3]);
    A4 = -einsum_kg(I4,t2b,'ncjk,abin->abcijk');
    A4 = A4 - permute(A4,[1,3,2,4,5,6]);
    A13 = A1 + A3;
    A13 = A13 - permute(A13,[1,3,2,4,5,6]) - permute(A13,[1,2,3,4,6,5]) + permute(A13,[1,3,2,4,6,5]);
    D12 = A13 + A2 + A4;

    %D13 = -einsum_kg(H3C.ovvovo,r2b,'mbciek,aemj->abcijk');
    %D13 = D13 - permute(D13,[1,2,3,4,6,5]);
    I1 = -einsum_kg(H2B.ovvv,r2b,'mbfe,aemj->abfj');
    A1 = einsum_kg(I1,t2b,'abfj,fcik->abcijk');
    A1 = A1 - permute(A1,[1,3,2,4,5,6]) - permute(A1,[1,2,3,4,6,5]) + permute(A1,[1,3,2,4,6,5]);
    I2 = -einsum_kg(H2B.ooov,r2b,'mnie,aemj->anij');
    A2 = -einsum_kg(I2,t2c,'anij,bcnk->abcijk');
    A2 = A2 - permute(A2,[1,2,3,4,6,5]);
    D13 = A1 + A2;

    %D14 = -einsum_kg(H3C.vvovoo,r2b,'acmekj,ebim->abcijk');
    %D14 = D14 - permute(D14,[1,3,2,4,5,6]);
    I1 = -einsum_kg(H2B.oovo,r2b,'nmej,ebim->nbij');
    A1 = -einsum_kg(I1,t2b,'nbij,acnk->abcijk');
    A1 = A1 - permute(A1,[1,3,2,4,5,6]) - permute(A1,[1,2,3,4,6,5]) + permute(A1,[1,3,2,4,6,5]);
    I2 = -einsum_kg(H2B.vovv,r2b,'amef,ebim->abif');
    A2 = einsum_kg(I2,t2c,'abif,fcjk->abcijk');
    A2 = A2 - permute(A2,[1,3,2,4,5,6]);
    D14 = A1 + A2;
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

    EOMMM23C =  D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 ...
              + D9 + D10 + D11 + D12 + D13 + D14 ...
              + D15 + D16 + D17 + D18 + D19 + D20 ...
              + D21 + D22 + D23;
     
    fprintf(' finished in %4.2f s\n',toc)
     
end