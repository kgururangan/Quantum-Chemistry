function [EOMMM23] = build_EOMMM23(t1,t2,r1,r2,HBar)

    H3 = HBar{3}; H2 = HBar{2}; 
    H_vvoooo = H3{2,2,1,1,1,1};
    H_vvooov = H3{2,2,1,1,1,2};
    H_vvvoov = H3{2,2,2,1,1,2};
    H_vvvvvo = H3{2,2,2,2,2,1};
    H_oovooo = H3{1,1,2,1,1,1};
    H_vooo = H2{2,1,1,1};
    H_vvov = H2{2,2,1,2};
    V = HBar{2}{1,1,2,2};
    
    fprintf('EOM-MMCC(2,3) construction... ')
    
    tic
    
    EOMMM23 = -einsum_kg(H_vvoooo,r1,'acmikj,bm->abcijk')...
              +einsum_kg(H_vvoooo,r1,'bcmikj,am->abcijk')... % (ab)
              +einsum_kg(H_vvoooo,r1,'abmikj,cm->abcijk');   % (bc)

    EOMMM23 = EOMMM23 +einsum_kg(H_vvvoov,r1,'acbike,ej->abcijk')...
                      -einsum_kg(H_vvvoov,r1,'acbije,ek->abcijk')... % (kj)
                      -einsum_kg(H_vvvoov,r1,'acbjke,ei->abcijk');   % (ij)
%     
    EOMMM23 = EOMMM23 +0.5*einsum_kg(H_vvvvvo,r2,'abcefk,efij->abcijk')...
                      -0.5*einsum_kg(H_vvvvvo,r2,'abcefi,efkj->abcijk')... % (ki)
                      -0.5*einsum_kg(H_vvvvvo,r2,'abcefj,efik->abcijk');   % (kj)
             
    EOMMM23 = EOMMM23 +0.5*einsum_kg(H_oovooo,r2,'mncijk,abmn->abcijk')... 
                      -0.5*einsum_kg(H_oovooo,r2,'mnaijk,cbmn->abcijk')... % (ac)
                      -0.5*einsum_kg(H_oovooo,r2,'mnbijk,acmn->abcijk');   % (bc)
%              
    EOMMM23 = EOMMM23 +einsum_kg(H_vvooov,r2,'abmije,ecmk->abcijk')... 
                      -einsum_kg(H_vvooov,r2,'cbmije,eamk->abcijk')... % (ac)
                      -einsum_kg(H_vvooov,r2,'acmije,ebmk->abcijk')... % (bc)
                      -einsum_kg(H_vvooov,r2,'abmkje,ecmi->abcijk')... % (ik)
                      -einsum_kg(H_vvooov,r2,'abmike,ecmj->abcijk')... % (jk)
                      +einsum_kg(H_vvooov,r2,'cbmkje,eami->abcijk')... % (ik)(ac)
                      +einsum_kg(H_vvooov,r2,'cbmike,eamj->abcijk')... % (jk)(ac)
                      +einsum_kg(H_vvooov,r2,'acmkje,ebmi->abcijk')... % (ik)(bc)
                      +einsum_kg(H_vvooov,r2,'acmike,ebmj->abcijk');   % (jk)(bc)
             
     EOMMM23 = EOMMM23 -einsum_kg(H_vooo,r2,'amij,bcmk->abcijk')...
                       +einsum_kg(H_vooo,r2,'amkj,bcmi->abcijk')... % (ki)
                       +einsum_kg(H_vooo,r2,'amik,bcmj->abcijk')... % (kj)
                       +einsum_kg(H_vooo,r2,'bmij,acmk->abcijk')... % (ab)
                       +einsum_kg(H_vooo,r2,'cmij,bamk->abcijk')... % (ac)
                       -einsum_kg(H_vooo,r2,'bmkj,acmi->abcijk')... % (ki)(ab)
                       -einsum_kg(H_vooo,r2,'cmkj,bami->abcijk')... % (ki)(ac)
                       -einsum_kg(H_vooo,r2,'bmik,acmj->abcijk')... % (kj)(ab)
                       -einsum_kg(H_vooo,r2,'cmik,bamj->abcijk');   % (kj)(ac)
             
     EOMMM23 = EOMMM23 +einsum_kg(H_vvov,r2,'abie,ecjk->abcijk')...
                       -einsum_kg(H_vvov,r2,'cbie,eajk->abcijk')... % (ac)
                       -einsum_kg(H_vvov,r2,'acie,ebjk->abcijk')... % (bc)
                       -einsum_kg(H_vvov,r2,'abke,ecji->abcijk')... % (ik)
                       -einsum_kg(H_vvov,r2,'abje,ecik->abcijk')... % (ij)
                       +einsum_kg(H_vvov,r2,'cbke,eaji->abcijk')... % (ac)(ik)
                       +einsum_kg(H_vvov,r2,'cbje,eaik->abcijk')... % (ac)(ij)
                       +einsum_kg(H_vvov,r2,'acke,ebji->abcijk')... % (bc)(ik)
                       +einsum_kg(H_vvov,r2,'acje,ebik->abcijk');   % (bc)(ij)

     I1 = einsum_kg(V,r1,'mnef,fn->me');
     EOMMM23 = EOMMM23 - einsum_kg(einsum_kg(I1,t2,'me,ecjk->mcjk'),t2,'mcjk,abim->abcijk') ...
                         + einsum_kg(einsum_kg(I1,t2,'me,ecik->mcik'),t2,'mcik,abjm->abcijk') ... % (ij) 
                         + einsum_kg(einsum_kg(I1,t2,'me,ecji->mcji'),t2,'mcji,abkm->abcijk') ... % (ik), only need A(i/jk) i think...
                         + einsum_kg(einsum_kg(I1,t2,'me,eajk->majk'),t2,'majk,cbim->abcijk') ... % (ac)
                         + einsum_kg(einsum_kg(I1,t2,'me,ebjk->mbjk'),t2,'mbjk,acim->abcijk') ... % (bc)
                         - einsum_kg(einsum_kg(I1,t2,'me,eaik->maik'),t2,'maik,cbjm->abcijk') ... % (ij)(ac)
                         - einsum_kg(einsum_kg(I1,t2,'me,ebik->mbik'),t2,'mbik,acjm->abcijk') ... % (ij)(bc)
                         - einsum_kg(einsum_kg(I1,t2,'me,eaji->maji'),t2,'maji,cbkm->abcijk') ... % (ik)(ac)
                         - einsum_kg(einsum_kg(I1,t2,'me,ebji->mbji'),t2,'mbji,ackm->abcijk');    % (ik)(bc)
%       
     
     fprintf(' finished in %4.2f s\n',toc)
             

end

