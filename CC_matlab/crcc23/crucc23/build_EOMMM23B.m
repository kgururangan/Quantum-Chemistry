function [EOMMM23B] = build_EOMMM23B(cc_t,HBar_t,iroot)

    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C; 
    %H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3D = HBar_t.H3D; H3C = HBar_t.H3C;
    
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;

    r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
    r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};
    
    fprintf('EOM-MMCC(2,3)B construction... ')
    
    tic
    
%     D1 = einsum_kg(H3B.vvvoov,r1b,'abcije,ek->abcijk');
    I1 = einsum_kg(H2B.ovov,r1b,'mcie,ek->mcik');
    I2 = einsum_kg(H2B.vvvv,r1b,'acfe,ek->acfk');
    I3 = einsum_kg(H2B.voov,r1b,'amie,ek->amik');
    A1 = -einsum_kg(I1,t2a,'mcik,abmj->abcijk');
    A1 = A1 - permute(A1,[1,2,3,5,4,6]);
    A2 = einsum_kg(I2,t2a,'acfk,fbij->abcijk');
    A2 = A2 - permute(A2,[2,1,3,4,5,6]);
    A3 = -einsum_kg(I3,t2b,'amik,bcjm->abcijk');
    A3 = A3 - permute(A3,[2,1,3,4,5,6]) - permute(A3,[1,2,3,5,4,6]) + permute(A3,[2,1,3,5,4,6]);
    D1 = A1 + A2 + A3;

%     D2 = einsum_kg(H3B.vvvovo,r1a,'abciek,ej->abcijk');
%     D2 = D2 - permute(D2,[1,2,3,5,4,6]);
    I1 = einsum_kg(H2B.ovvo,r1a,'mcek,ej->mcjk');
    I2 = einsum_kg(H2B.vovo,r1a,'bmek,ej->bmjk');
    I3 = einsum_kg(H2A.voov,r1a,'amie,ej->amij'); I3 = I3 - permute(I3,[1,2,4,3]);
    I4 = einsum_kg(H2B.vvvv,r1a,'bcef,ej->bcjf');
    I5 = einsum_kg(H2A.vvvv,r1a,'abfe,ej->abfj');
    A1 = -einsum_kg(I1,t2a,'mcjk,abim->abcijk');
    A2 = -einsum_kg(I2,t2b,'bmjk,acim->abcijk');
    A3 = -einsum_kg(I3,t2b,'amij,bcmk->abcijk');
    A4 = einsum_kg(I4,t2b,'bcjf,afik->abcijk');
    A5 = einsum_kg(I5,t2b,'abfj,fcik->abcijk');
    B1 = A1 + A5; 
    B2 = A2 + A4;
    B1 = B1 - permute(B1,[1,2,3,5,4,6]);
    B2 = B2 - permute(B2,[1,2,3,5,4,6]) - permute(B2,[2,1,3,4,5,6]) + permute(B2,[2,1,3,5,4,6]);
    A3 = A3 - permute(A3,[2,1,3,4,5,6]);
    D2 = B1 + B2 + A3;
    
%     D3 = -einsum_kg(H3B.vvoooo,r1b,'abmijk,cm->abcijk');
    I1 = -einsum_kg(H2B.vovo,r1b,'bmek,cm->bcek');
    I2 = -einsum_kg(H2B.voov,r1b,'bmje,cm->bcje');
    I3 = -einsum_kg(H2B.oooo,r1b,'nmjk,cm->ncjk');
    A1 = einsum_kg(I1,t2a,'bcek,aeij->abcijk');
    A2 = einsum_kg(I2,t2b,'bcje,aeik->abcijk');
    A3 = -einsum_kg(I3,t2a,'ncjk,abin->abcijk');
    A1 = A1 - permute(A1,[2,1,3,4,5,6]);
    A2 = A2 - permute(A2,[2,1,3,4,5,6]) - permute(A2,[1,2,3,5,4,6]) + permute(A2,[2,1,3,5,4,6]);
    A3 = A3 - permute(A3,[1,2,3,5,4,6]);
    D3 = A1 + A2 + A3;

%     D4 = -einsum_kg(H3B.vovooo,r1a,'amcijk,bm->abcijk');
%     D4 = D4 - permute(D4,[2,1,3,4,5,6]);
    I1 = -einsum_kg(H2B.ovov,r1a,'mcje,bm->bcje');
    I2 = -einsum_kg(H2B.ovvo,r1a,'mcek,bm->bcek');
    I3 = -einsum_kg(H2B.oooo,r1a,'mnjk,bm->bnjk');
    I4 = -einsum_kg(H2A.oooo,r1a,'mnji,bm->bnji');
    I5 = einsum_kg(H2A.voov,r1a,'amje,bm->abej'); I5 = I5 - permute(I5,[2,1,3,4]);
    A1 = einsum_kg(I1,t2b,'bcje,aeik->abcijk');
    A2 = einsum_kg(I2,t2a,'bcek,aeij->abcijk');
    A3 = -einsum_kg(I3,t2b,'bnjk,acin->abcijk');
    A4 = -einsum_kg(I4,t2b,'bnji,acnk->abcijk');
    A5 = einsum_kg(I5,t2b,'abej,ecik->abcijk');
    B1 = A1 + A3;
    B1 = B1 - permute(B1,[2,1,3,4,5,6]) - permute(B1,[1,2,3,5,4,6]) + permute(B1,[2,1,3,5,4,6]);
    B2 = A2 + A4;
    B2 = B2 - permute(B2,[2,1,3,4,5,6]);
    A5 = A5 - permute(A5,[1,2,3,5,4,6]);
    D4 = B1 + B2 + A5;

%     D5 = 0.5*einsum_kg(H3B.oovooo,r2a,'mncijk,abmn->abcijk');
    I1 = 0.5*einsum_kg(H2A.oovo,r2a,'mnej,abmn->abej');
    D5 = einsum_kg(I1,t2b,'abej,ecik->abcijk');
    D5 = D5 - permute(D5,[1,2,3,5,4,6]);

%     D6 = 0.5*einsum_kg(H3B.vvvvvo,r2a,'abcefk,efij->abcijk');
    I1 = 0.5*einsum_kg(H2A.vovv,r2a,'bmef,efji->mbij');
    D6 = -einsum_kg(I1,t2b,'mbij,acmk->abcijk');
    D6 = D6 - permute(D6,[2,1,3,4,5,6]);

%     D7 = einsum_kg(H3B.vooooo,r2b,'bmnjik,acmn->abcijk');
%     D7 = D7 - permute(D7,[2,1,3,4,5,6]);
    I1 = einsum_kg(H2B.oovo,r2b,'mnek,acmn->acek');
    I2 = einsum_kg(H2B.ooov,r2b,'mnie,acmn->acie');
    A1 = einsum_kg(I1,t2a,'acek,ebij->abcijk');
    A1 = A1 - permute(A1,[2,1,3,4,5,6]);
    A2 = einsum_kg(I2,t2b,'acie,bejk->abcijk');
    A2 = A2 - permute(A2,[2,1,3,4,5,6]) - permute(A2,[1,2,3,5,4,6]) + permute(A2,[2,1,3,5,4,6]);
    D7 = A1 + A2;
    
%     D8 = einsum_kg(H3B.vvvovv,r2b,'bacjef,efik->abcijk');
%     D8 = D8 - permute(D8,[1,2,3,5,4,6]);
    I1 = einsum_kg(H2B.ovvv,r2b,'mcef,efik->mcik');
    I2 = einsum_kg(H2B.vovv,r2b,'amef,efik->amik');
    A1 = -einsum_kg(I1,t2a,'mcik,abmj->abcijk');
    A1 = A1 - permute(A1,[1,2,3,5,4,6]);
    A2 = -einsum_kg(I2,t2b,'amik,bcjm->abcijk');
    A2 = A2 - permute(A2,[2,1,3,4,5,6]) - permute(A2,[1,2,3,5,4,6]) + permute(A2,[2,1,3,5,4,6]);
    D8 = A1 + A2;

%     D9 = einsum_kg(H3A.vvooov,r2b,'abmije,ecmk->abcijk');
    I1 = einsum_kg(H2A.ooov,r2b,'nmie,ecmk->ncik');
    I2 = einsum_kg(H2A.vovv,r2b,'amfe,ecmk->acfk');
    A1 = -einsum_kg(I1,t2a,'ncik,abnj->abcijk');
    A2 = einsum_kg(I2,t2a,'acfk,fbij->abcijk');
    A1 = A1 - permute(A1,[1,2,3,5,4,6]);
    A2 = A2 - permute(A2,[2,1,3,4,5,6]);
    D9 = A1 + A2;

%     D10 = einsum_kg(H3B.vvooov,r2c,'abmije,ecmk->abcijk');
    I1 = einsum_kg(H2B.ooov,r2c,'nmie,ecmk->ncik');
    I2 = einsum_kg(H2B.vovv,r2c,'amfe,ecmk->acfk');
    A1 = -einsum_kg(I1,t2a,'ncik,abnj->abcijk');
    A2 = einsum_kg(I2,t2a,'acfk,fbij->abcijk');
    A1 = A1 - permute(A1,[1,2,3,5,4,6]);
    A2 = A2 - permute(A2,[2,1,3,4,5,6]);
    D10 = A1 + A2;

%     D11 = einsum_kg(H3B.vovovo,r2a,'amciek,bejm->abcijk');
%     D11 = D11 - permute(D11,[2,1,3,4,5,6]) - permute(D11,[1,2,3,5,4,6]) + permute(D11,[2,1,3,5,4,6]);
    I1 = einsum_kg(H2A.ooov,r2a,'nmie,ebmj->nbij'); I1 = I1 - permute(I1,[1,2,4,3]);
    I2 = einsum_kg(H2B.oovo,r2a,'mnek,bejm->bnjk');
    I3 = einsum_kg(H2A.vovv,r2a,'amfe,bejm->abfj'); I3 = I3 - permute(I3,[2,1,3,4]);
    I4 = einsum_kg(H2B.ovvv,r2a,'mcef,bejm->bcjf');
    A1 = -einsum_kg(I1,t2b,'nbij,acnk->abcijk');
    A1 = A1 - permute(A1,[2,1,3,4,5,6]);
    A2 = -einsum_kg(I2,t2b,'bnjk,acin->abcijk');
    A3 = einsum_kg(I3,t2b,'abfj,fcik->abcijk');
    A3 = A3 - permute(A3,[1,2,3,5,4,6]);
    A4 = einsum_kg(I4,t2b,'bcjf,afik->abcijk');
    B1 = A2 + A4;
    B1 = B1 - permute(B1,[2,1,3,4,5,6]) - permute(B1,[1,2,3,5,4,6]) + permute(B1,[2,1,3,5,4,6]);
    D11 = A1 + A3 + B1;

%     D12 = einsum_kg(H3C.vovovo,r2b,'amciek,bejm->abcijk');
%     D12 = D12 - permute(D12,[2,1,3,4,5,6]) - permute(D12,[1,2,3,5,4,6]) + permute(D12,[2,1,3,5,4,6]);
    I1 = einsum_kg(H2B.ooov,r2b,'nmie,bejm->nbij'); I1 = I1 - permute(I1,[1,2,4,3]);
    I2 = einsum_kg(H2B.vovv,r2b,'amfe,bejm->abfj'); I2 = I2 - permute(I2,[2,1,3,4]);
    I3 = einsum_kg(H2C.vovv,r2b,'cmfe,bejm->bcjf');
    I4 = einsum_kg(H2C.ooov,r2b,'nmke,bejm->bnjk');
    A1 = -einsum_kg(I1,t2b,'nbij,acnk->abcijk');
    A1 = A1 - permute(A1,[2,1,3,4,5,6]);
    A2 = einsum_kg(I2,t2b,'abfj,fcik->abcijk');
    A2 = A2 - permute(A2,[1,2,3,5,4,6]);
    A3 = einsum_kg(I3,t2b,'bcjf,afik->abcijk');
    A4 = -einsum_kg(I4,t2b,'bnjk,acin->abcijk');
    B1 = A3 + A4;
    B1 = B1 - permute(B1,[2,1,3,4,5,6]) - permute(B1,[1,2,3,5,4,6]) + permute(B1,[2,1,3,5,4,6]);
    D12 = A1 + A2 + B1;

%     D13 = -einsum_kg(H3B.vvoovo,r2b,'banjek,ecin->abcijk');
%     D13 = D13 - permute(D13,[1,2,3,5,4,6]);
    I1 = -einsum_kg(H2B.oovo,r2b,'nmek,ecim->ncik');
    I2 = -einsum_kg(H2B.vovv,r2b,'amef,ecim->acif');
    A1 = -einsum_kg(I1,t2a,'ncik,abnj->abcijk');
    A1 = A1 - permute(A1,[1,2,3,5,4,6]);
    A2 = einsum_kg(I2,t2b,'acif,bfjk->abcijk');
    A2 = A2 - permute(A2,[2,1,3,4,5,6]) - permute(A2,[1,2,3,5,4,6]) + permute(A2,[2,1,3,5,4,6]);
    D13 = A1 + A2;

%     D14 = -einsum_kg(H3B.vovoov,r2b,'bmcjie,aemk->abcijk');
%     D14 = D14 - permute(D14,[2,1,3,4,5,6]);
    I1 = -einsum_kg(H2B.ooov,r2b,'mnie,aemk->anik');
    I2 = -einsum_kg(H2B.ovvv,r2b,'mcfe,aemk->acfk');
    A1 = -einsum_kg(I1,t2b,'anik,bcjn->abcijk');
    A1 = A1 - permute(A1,[2,1,3,4,5,6]) - permute(A1,[1,2,3,5,4,6]) + permute(A1,[2,1,3,5,4,6]);
    A2 = einsum_kg(I2,t2a,'acfk,fbij->abcijk');
    A2 = A2 - permute(A2,[2,1,3,4,5,6]);
    D14 = A1 + A2;

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

    I1A = einsum_kg(H2A.oovv,r1a,'mnef,fn->me') + einsum_kg(H2B.oovv,r1b,'mnef,fn->me');
    I1B = einsum_kg(H2C.oovv,r1b,'nmfe,fn->me') + einsum_kg(H2B.oovv,r1a,'nmfe,fn->me');

    D21 = -einsum_kg(einsum_kg(I1A,t2a,'me,aeij->amij'),t2b,'amij,bcmk->abcijk');
    D21 = D21 - permute(D21,[2,1,3,4,5,6]);

    D22 = -einsum_kg(einsum_kg(I1A,t2a,'me,abim->abie'),t2b,'abie,ecjk->abcijk');
    D22 = D22 - permute(D22,[1,2,3,5,4,6]);

    D23 = -einsum_kg(einsum_kg(I1B,t2b,'me,aeik->amik'),t2b,'amik,bcjm->abcijk');
    D23 = D23 - permute(D23,[2,1,3,4,5,6]) - permute(D23,[1,2,3,5,4,6]) + permute(D23,[2,1,3,5,4,6]);

    EOMMM23B =  D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 ...
              + D9 + D10 + D11 + D12 + D13 + D14 ...
              + D15 + D16 + D17 + D18 + D19 + D20 ...
              + D21 + D22 + D23;
     
    fprintf(' finished in %4.2f s\n',toc)
     
end