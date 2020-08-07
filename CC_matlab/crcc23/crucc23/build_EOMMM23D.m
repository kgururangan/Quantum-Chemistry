function [EOMMM23D] = build_EOMMM23D(cc_t,HBar_t,iroot)

    H3D = HBar_t.H3D; H2C = HBar_t.H2C; H3C = HBar_t.H3C;
    VC = HBar_t.H2C.oovv; VB = HBar_t.H2B.oovv;
    t2c = cc_t.t2c;

    r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
    r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};
    
    fprintf('EOM-MMCC(2,3)D construction... ')
    
    tic
    
    EOMMM23D = -einsum_kg(H3D.vvoooo,r1b,'acmikj,bm->abcijk')...
               +einsum_kg(H3D.vvoooo,r1b,'bcmikj,am->abcijk')... % (ab)
               +einsum_kg(H3D.vvoooo,r1b,'abmikj,cm->abcijk');   % (bc)
%                 
    EOMMM23D = EOMMM23D +einsum_kg(H3D.vvvoov,r1b,'acbike,ej->abcijk')...
                        -einsum_kg(H3D.vvvoov,r1b,'acbije,ek->abcijk')... % (kj)
                        -einsum_kg(H3D.vvvoov,r1b,'acbjke,ei->abcijk');   % (ij)
%               
    EOMMM23D = EOMMM23D +0.5*einsum_kg(H3D.vvvvvo,r2c,'abcefk,efij->abcijk')...
                        -0.5*einsum_kg(H3D.vvvvvo,r2c,'abcefi,efkj->abcijk')... % (ki)
                        -0.5*einsum_kg(H3D.vvvvvo,r2c,'abcefj,efik->abcijk');   % (kj)
%             
    EOMMM23D = EOMMM23D +0.5*einsum_kg(H3D.oovooo,r2c,'mncijk,abmn->abcijk')... 
                        -0.5*einsum_kg(H3D.oovooo,r2c,'mnaijk,cbmn->abcijk')... % (ac)
                        -0.5*einsum_kg(H3D.oovooo,r2c,'mnbijk,acmn->abcijk');   % (bc)
%             
    EOMMM23D = EOMMM23D +einsum_kg(H3D.vvooov,r2c,'abmije,ecmk->abcijk')... 
                        -einsum_kg(H3D.vvooov,r2c,'cbmije,eamk->abcijk')... % (ac)
                        -einsum_kg(H3D.vvooov,r2c,'acmije,ebmk->abcijk')... % (bc)
                        -einsum_kg(H3D.vvooov,r2c,'abmkje,ecmi->abcijk')... % (ik)
                        -einsum_kg(H3D.vvooov,r2c,'abmike,ecmj->abcijk')... % (jk)
                        +einsum_kg(H3D.vvooov,r2c,'cbmkje,eami->abcijk')... % (ik)(ac)
                        +einsum_kg(H3D.vvooov,r2c,'cbmike,eamj->abcijk')... % (jk)(ac)
                        +einsum_kg(H3D.vvooov,r2c,'acmkje,ebmi->abcijk')... % (ik)(bc)
                        +einsum_kg(H3D.vvooov,r2c,'acmike,ebmj->abcijk');   % (jk)(bc)
                    
     % PROBLEM - should adapt indices to D case properly!               
     EOMMM23D = EOMMM23D +einsum_kg(H3C.ovvvoo,r2b,'mabeij,ecmk->abcijk')... 
                         -einsum_kg(H3C.ovvvoo,r2b,'mcbeij,eamk->abcijk')... % (ac)
                         -einsum_kg(H3C.ovvvoo,r2b,'maceij,ebmk->abcijk')... % (bc)
                         -einsum_kg(H3C.ovvvoo,r2b,'mabekj,ecmi->abcijk')... % (ik)
                         -einsum_kg(H3C.ovvvoo,r2b,'mabeik,ecmj->abcijk')... % (jk)
                         +einsum_kg(H3C.ovvvoo,r2b,'mcbekj,eami->abcijk')... % (ik)(ac)
                         +einsum_kg(H3C.ovvvoo,r2b,'mcbeik,eamj->abcijk')... % (jk)(ac)
                         +einsum_kg(H3C.ovvvoo,r2b,'macekj,ebmi->abcijk')... % (ik)(bc)
                         +einsum_kg(H3C.ovvvoo,r2b,'maceik,ebmj->abcijk');   % (jk)(bc)
            
     EOMMM23D = EOMMM23D -einsum_kg(H2C.vooo,r2c,'amij,bcmk->abcijk')...
                         +einsum_kg(H2C.vooo,r2c,'amkj,bcmi->abcijk')... % (ki)
                         +einsum_kg(H2C.vooo,r2c,'amik,bcmj->abcijk')... % (kj)
                         +einsum_kg(H2C.vooo,r2c,'bmij,acmk->abcijk')... % (ab)
                         +einsum_kg(H2C.vooo,r2c,'cmij,bamk->abcijk')... % (ac)
                         -einsum_kg(H2C.vooo,r2c,'bmkj,acmi->abcijk')... % (ki)(ab)
                         -einsum_kg(H2C.vooo,r2c,'cmkj,bami->abcijk')... % (ki)(ac)
                         -einsum_kg(H2C.vooo,r2c,'bmik,acmj->abcijk')... % (kj)(ab)
                         -einsum_kg(H2C.vooo,r2c,'cmik,bamj->abcijk');   % (ij)(ac)
             
     EOMMM23D = EOMMM23D +einsum_kg(H2C.vvov,r2c,'abie,ecjk->abcijk')...
                         -einsum_kg(H2C.vvov,r2c,'cbie,eajk->abcijk')... % (ac)
                         -einsum_kg(H2C.vvov,r2c,'acie,ebjk->abcijk')... % (bc)
                         -einsum_kg(H2C.vvov,r2c,'abke,ecji->abcijk')... % (ik)
                         -einsum_kg(H2C.vvov,r2c,'abje,ecik->abcijk')... % (ij)
                         +einsum_kg(H2C.vvov,r2c,'cbke,eaji->abcijk')... % (ac)(ik)
                         +einsum_kg(H2C.vvov,r2c,'cbje,eaik->abcijk')... % (ac)(ij)
                         +einsum_kg(H2C.vvov,r2c,'acke,ebji->abcijk')... % (bc)(ik)
                         +einsum_kg(H2C.vvov,r2c,'acje,ebik->abcijk');   % (bc)(ij)
                     
     I1 = einsum_kg(VB,r1a,'mnef,fn->me');
     EOMMM23D = EOMMM23D - einsum_kg(einsum_kg(I1,t2c,'me,ecjk->mcjk'),t2c,'mcjk,abim->abcijk') ...
                         + einsum_kg(einsum_kg(I1,t2c,'me,ecik->mcik'),t2c,'mcik,abjm->abcijk') ... % (ij) 
                         + einsum_kg(einsum_kg(I1,t2c,'me,ecji->mcji'),t2c,'mcji,abkm->abcijk') ... % (ik), only need A(i/jk) i think...
                         + einsum_kg(einsum_kg(I1,t2c,'me,eajk->majk'),t2c,'majk,cbim->abcijk') ... % (ac)
                         + einsum_kg(einsum_kg(I1,t2c,'me,ebjk->mbjk'),t2c,'mbjk,acim->abcijk') ... % (bc)
                         - einsum_kg(einsum_kg(I1,t2c,'me,eaik->maik'),t2c,'maik,cbjm->abcijk') ... % (ij)(ac)
                         - einsum_kg(einsum_kg(I1,t2c,'me,ebik->mbik'),t2c,'mbik,acjm->abcijk') ... % (ij)(bc)
                         - einsum_kg(einsum_kg(I1,t2c,'me,eaji->maji'),t2c,'maji,cbkm->abcijk') ... % (ik)(ac)
                         - einsum_kg(einsum_kg(I1,t2c,'me,ebji->mbji'),t2c,'mbji,ackm->abcijk');    % (ik)(bc)
     
     I2 = einsum_kg(VC,r1b,'mnef,fn->me');
     EOMMM23D = EOMMM23D - einsum_kg(einsum_kg(I2,t2c,'me,ecjk->mcjk'),t2c,'mcjk,abim->abcijk') ...
                         + einsum_kg(einsum_kg(I2,t2c,'me,ecik->mcik'),t2c,'mcik,abjm->abcijk') ... % (ij) 
                         + einsum_kg(einsum_kg(I2,t2c,'me,ecji->mcji'),t2c,'mcji,abkm->abcijk') ... % (ik), only need A(i/jk) i think...
                         + einsum_kg(einsum_kg(I2,t2c,'me,eajk->majk'),t2c,'majk,cbim->abcijk') ... % (ac)
                         + einsum_kg(einsum_kg(I2,t2c,'me,ebjk->mbjk'),t2c,'mbjk,acim->abcijk') ... % (bc)
                         - einsum_kg(einsum_kg(I2,t2c,'me,eaik->maik'),t2c,'maik,cbjm->abcijk') ... % (ij)(ac)
                         - einsum_kg(einsum_kg(I2,t2c,'me,ebik->mbik'),t2c,'mbik,acjm->abcijk') ... % (ij)(bc)
                         - einsum_kg(einsum_kg(I2,t2c,'me,eaji->maji'),t2c,'maji,cbkm->abcijk') ... % (ik)(ac)
                         - einsum_kg(einsum_kg(I2,t2c,'me,ebji->mbji'),t2c,'mbji,ackm->abcijk');    % (ik)(bc)
     
%      I1 = einsum_kg(VB,r1a,'mnef,fn->me');
%      EOMMM23D = EOMMM23D - einsum_kg(einsum_kg(I1,t2c,'me,ecjk->mcjk'),t2c,'mcjk,abim->abcijk') ...
%                          + einsum_kg(einsum_kg(I1,t2c,'me,ecik->mcik'),t2c,'mcik,abjm->abcijk') ... % (ij) 
%                          + einsum_kg(einsum_kg(I1,t2c,'me,ecji->mcji'),t2c,'mcji,abkm->abcijk');    % (ik) 
%      
%      I2 = einsum_kg(VC,r1b,'mnef,fn->me');
%      EOMMM23D = EOMMM23D - einsum_kg(einsum_kg(I2,t2c,'me,ecjk->mcjk'),t2c,'mcjk,abim->abcijk') ...
%                          + einsum_kg(einsum_kg(I2,t2c,'me,ecik->mcik'),t2c,'mcik,abjm->abcijk') ... % (ij)
%                          + einsum_kg(einsum_kg(I2,t2c,'me,ecji->mcji'),t2c,'mcji,abkm->abcijk');    % (ik)
     
     %EOMMM23D = EOMMM23D + r0*MM23D
     
     fprintf(' finished in %4.2f s\n',toc)
     
end