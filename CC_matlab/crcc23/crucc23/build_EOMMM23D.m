function [EOMMM23D] = build_EOMMM23D(cc_t,HBar_t,iroot)

    H3D = HBar_t.H3D; H2C = HBar_t.H2C; H3C = HBar_t.H3C;
    VC = HBar_t.H2C.oovv; VB = HBar_t.H2B.oovv;
    t2c = cc_t.t2c;

    r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
    r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};
    
    fprintf('EOM-MMCC(2,3)D construction... ')
    
    tic
    
    D1 = -einsum_kg(H3D.vvoooo,r1b,'acmikj,bm->abcijk');
    D1 = D1 - permute(D1,[2,1,3,4,5,6]) - permute(D1,[1,3,2,4,5,6]);
    
    D2 = +einsum_kg(H3D.vvvoov,r1b,'acbike,ej->abcijk');
    D2 = D2 - permute(D2,[1,2,3,4,6,5]) - permute(D2,[1,2,3,5,4,6]);
    
    D3 = +0.5*einsum_kg(H3D.vvvvvo,r2c,'abcefk,efij->abcijk');
    D3 = D3 - permute(D3,[1,2,3,6,5,4]) - permute(D3,[1,2,3,4,6,5]);
    
    D4 = +0.5*einsum_kg(H3D.oovooo,r2c,'mncijk,abmn->abcijk');
    D4 = D4 - permute(D4,[3,2,1,4,5,6]) - permute(D4,[1,3,2,4,5,6]);
    
    D5 = +einsum_kg(H3D.vvooov,r2c,'abmije,ecmk->abcijk');
    D5 = D5 - permute(D5,[3,2,1,4,5,6]) - permute(D5,[1,3,2,4,5,6]);
    D5 = D5 - permute(D5,[1,2,3,4,6,5]) - permute(D5,[1,2,3,6,5,4]);
    
    D6 = +einsum_kg(H3C.ovvvoo,r2b,'mabeij,ecmk->abcijk');
    D6 = D6 - permute(D6,[3,2,1,4,5,6]) - permute(D6,[1,3,2,4,5,6]);
    D6 = D6 - permute(D6,[1,2,3,6,5,4]) - permute(D6,[1,2,3,4,6,5]);
    
    D7 = -einsum_kg(H2C.vooo,r2c,'amij,bcmk->abcijk');
    D7 = D7 - permute(D7,[2,1,3,4,5,6]) - permute(D7,[3,2,1,4,5,6]);
    D7 = D7 - permute(D7,[1,2,3,6,5,4]) - permute(D7,[1,2,3,4,6,5]);
    
    D8 = +einsum_kg(H2C.vvov,r2c,'abie,ecjk->abcijk');
    D8 = D8 - permute(D8,[3,2,1,4,5,6]) - permute(D8,[1,3,2,4,5,6]);
    D8 = D8 - permute(D8,[1,2,3,6,5,4]) - permute(D8,[1,2,3,5,4,6]);
    
    I1 = einsum_kg(VB,r1a,'mnef,fn->me');
    D9 = -einsum_kg(einsum_kg(I1,t2c,'me,ecjk->mcjk'),t2c,'mcjk,abim->abcijk');
    D9 = D9 - permute(D9,[3,2,1,4,5,6]) - permute(D9,[1,3,2,4,5,6]);
    D9 = D9 - permute(D9,[1,2,3,5,4,6]) - permute(D9,[1,2,3,6,5,4]);
    
    I2 = einsum_kg(VC,r1b,'mnef,fn->me');
    D10 =  -einsum_kg(einsum_kg(I2,t2c,'me,ecjk->mcjk'),t2c,'mcjk,abim->abcijk');
    D10 = D10 - permute(D10,[3,2,1,4,5,6]) - permute(D10,[1,3,2,4,5,6]);
    D10 = D10 - permute(D10,[1,2,3,5,4,6]) - permute(D10,[1,2,3,6,5,4]);
    
    EOMMM23D = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10;
     
    fprintf(' finished in %4.2f s\n',toc)
     
end