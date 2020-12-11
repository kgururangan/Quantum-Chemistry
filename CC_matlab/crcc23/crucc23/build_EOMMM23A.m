function [EOMMM23A] = build_EOMMM23A(cc_t,HBar_t,iroot)

    H2A = HBar_t.H2A; H2B = HBar_t.H2B; 

    t2a = cc_t.t2a;
    r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
    r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; 
    
    fprintf('\nEOM-MMCC(2,3)A construction... ')
    
    tic
    
%     D1 = -einsum_kg(H3A.vvoooo,r1a,'acmikj,bm->abcijk');
%     D1 = D1 - permute(D1,[2,1,3,4,5,6]) - permute(D1,[1,3,2,4,5,6]);
    
    I1 = -einsum_kg(H2A.voov,r1a,'amie,cm->acie');
    I1 = I1 - permute(I1,[2,1,3,4]);
    I2 = -einsum_kg(H2A.oooo,r1a,'nmjk,cm->ncjk');
    D1 = einsum_kg(I1,t2a,'abie,ecjk->abcijk') - einsum_kg(I2,t2a,'ncjk,abin->abcijk');
    D1 = D1 - permute(D1,[3,2,1,4,5,6]) - permute(D1,[1,3,2,4,5,6]);
    D1 = D1 - permute(D1,[1,2,3,6,5,4]) - permute(D1,[1,2,3,5,4,6]);
    
%     D2 = +einsum_kg(H3A.vvvoov,r1a,'acbike,ej->abcijk');
%     D2 = D2 - permute(D2,[1,2,3,4,6,5]) - permute(D2,[1,2,3,5,4,6]);

    I1 = einsum_kg(H2A.voov,r1a,'amie,ej->amij');
    I1 = I1 - permute(I1,[1,2,4,3]);
    I2 = einsum_kg(H2A.vvvv,r1a,'abfe,fi->abie');
    D2 = -einsum_kg(I1,t2a,'amij,bcmk->abcijk') + einsum_kg(I2,t2a,'cbke,aeij->abcijk');
    D2 = D2 - permute(D2,[2,1,3,4,5,6]) - permute(D2,[3,2,1,4,5,6]);
    D2 = D2 - permute(D2,[1,2,3,6,5,4]) - permute(D2,[1,2,3,4,6,5]);
    
%     D3 = +0.5*einsum_kg(H3A.vvvvvo,r2a,'abcefk,efij->abcijk');
%     D3 = D3 - permute(D3,[1,2,3,6,5,4]) - permute(D3,[1,2,3,4,6,5]);

    I1 = 0.5*einsum_kg(H2A.vovv,r2a,'amef,efij->amij');
    D3 = -einsum_kg(I1,t2a,'amij,bcmk->abcijk');
    D3 = D3 - permute(D3,[2,1,3,4,5,6]) - permute(D3,[3,2,1,4,5,6]);
    D3 = D3 - permute(D3,[1,2,3,6,5,4]) - permute(D3,[1,2,3,4,6,5]);
    
%     D4 = +0.5*einsum_kg(H3A.oovooo,r2a,'mncijk,abmn->abcijk');
%     D4 = D4 - permute(D4,[3,2,1,4,5,6]) - permute(D4,[1,3,2,4,5,6]);

    I1 = 0.5*einsum_kg(H2A.ooov,r2a,'mnie,abmn->abie');
    D4 = einsum_kg(I1,t2a,'abie,ecjk->abcijk');
    D4 = D4 - permute(D4,[3,2,1,4,5,6]) - permute(D4,[1,3,2,4,5,6]);
    D4 = D4 - permute(D4,[1,2,3,5,4,6]) - permute(D4,[1,2,3,6,5,4]);
    
%     D5 = +einsum_kg(H3A.vvooov,r2a,'abmije,ecmk->abcijk');
%     D5 = D5 - permute(D5,[3,2,1,4,5,6]) - permute(D5,[1,3,2,4,5,6]);
%     D5 = D5 - permute(D5,[1,2,3,6,5,4]) - permute(D5,[1,2,3,4,6,5]);

    I1 = einsum_kg(H2A.vovv,r2a,'bmfe,aeim->abif')...
            - einsum_kg(H2A.vovv,r2a,'amfe,beim->abif');
    I2 = einsum_kg(H2A.ooov,r2a,'nmje,cekm->cnkj')...
            - einsum_kg(H2A.ooov,r2a,'nmke,cejm->cnkj');
    D5 = einsum_kg(I1,t2a,'abif,fcjk->abcijk') - einsum_kg(I2,t2a,'cnkj,abin->abcijk');
    D5 = D5 - permute(D5,[3,2,1,4,5,6]) - permute(D5,[1,3,2,4,5,6]);
    D5 = D5 - permute(D5,[1,2,3,5,4,6]) - permute(D5,[1,2,3,6,5,4]);
    
%     D6 = +einsum_kg(H3B.vvooov,r2b,'abmije,cekm->abcijk');
%     D6 = D6 - permute(D6,[3,2,1,4,5,6]) - permute(D6,[1,3,2,4,5,6]);
%     D6 = D6 - permute(D6,[1,2,3,6,5,4]) - permute(D6,[1,2,3,4,6,5]);
    
    I1 = einsum_kg(H2B.vovv,r2b,'bmfe,aeim->abif')...
            - einsum_kg(H2B.vovv,r2b,'amfe,beim->abif');
    I2 = einsum_kg(H2B.ooov,r2b,'nmje,cekm->cnkj')...
            - einsum_kg(H2B.ooov,r2b,'nmke,cejm->cnkj');
    D6 = einsum_kg(I1,t2a,'abif,fcjk->abcijk') - einsum_kg(I2,t2a,'cnkj,abin->abcijk');
    D6 = D6 - permute(D6,[3,2,1,4,5,6]) - permute(D6,[1,3,2,4,5,6]);
    D6 = D6 - permute(D6,[1,2,3,5,4,6]) - permute(D6,[1,2,3,6,5,4]);
    
    D7 = -einsum_kg(H2A.vooo,r2a,'amij,bcmk->abcijk');
    D7 = D7 - permute(D7,[1,2,3,6,5,4]) - permute(D7,[1,2,3,4,6,5]);
    D7 = D7 - permute(D7,[2,1,3,4,5,6]) - permute(D7,[3,2,1,4,5,6]);
    
    D8 = +einsum_kg(H2A.vvov,r2a,'abie,ecjk->abcijk');
    D8 = D8 - permute(D8,[3,2,1,4,5,6]) - permute(D8,[1,3,2,4,5,6]);
    D8 = D8 - permute(D8,[1,2,3,5,4,6]) - permute(D8,[1,2,3,6,5,4]);
    
    I1 = einsum_kg(H2A.oovv,r1a,'mnef,fn->me') + einsum_kg(H2B.oovv,r1b,'mnef,fn->me');
    
    D9 = -einsum_kg(einsum_kg(I1,t2a,'me,ecjk->mcjk'),t2a,'mcjk,abim->abcijk');
    D9 = D9 - permute(D9,[1,2,3,5,4,6]) - permute(D9,[1,2,3,6,5,4]);
    D9 = D9 - permute(D9,[1,3,2,4,5,6]) - permute(D9,[3,2,1,4,5,6]);
    
    EOMMM23A = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9;
     
    fprintf(' finished in %4.2f s\n',toc)
             
end


