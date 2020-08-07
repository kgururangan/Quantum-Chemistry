function [EOMMM23A] = build_EOMMM23A(cc_t,HBar_t,iroot)

    H3A = HBar_t.H3A; H2A = HBar_t.H2A; H3B = HBar_t.H3B;
    VA = HBar_t.H2A.oovv; VB = HBar_t.H2B.oovv;
    t2a = cc_t.t2a;
    
    r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
    r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};
    
    fprintf('\nEOM-MMCC(2,3)A construction... ')
    
    tic
    
%      I_ov = -einsum_kg(H2A.oovv,r1a,'miae,em->ia') - einsum_kg(H2B.oovv,r1b,'ijab,bj->ia');
%      
%      I_ovoo = einsum_kg(I_ov,t2a,'ie,aekj->iajk')...
%               +0.5*einsum_kg(H2A.vovv,r2a,'aief,efkj->iajk')...
%               -einsum_kg(H2A.oooo,r1a,'imjk,am->iajk')...
%               +einsum_kg(H2A.ooov,r2a,'imje,aekm->iajk')-einsum_kg(H2A.voov,r1a,'aije,ek->iajk')...
%               -einsum_kg(H2A.ooov,r2a,'imke,aejm->iajk')+einsum_kg(H2A.voov,r1a,'aike,ej->iajk');
%      
%      I_vvov = 0.5*einsum_kg(H2A.vvvv,r1a,'abec,ei->abic')...
%               -einsum_kg(H2A.vovv,r2a,'bmec,aeim->abic')...
%               +0.25*einsum_kg(H2A.ooov,r2,'mnic,abmn->abic')...
%               +einsum_kg(H2A.voov,r1,'bmic,am->abic');
%           
%      EOMMM23A = 0.5*(einsum_kg(H2A.vvov,r2a,'abie,ecjk->abcijk')...
%                     -einsum_kg(H2A.vvov,r2a,'abje,ecik->abcijk')...
%                     -einsum_kg(H2A.vvov,r2a,'abke,ecji->abcijk'));
%      
%      EOMMM23A = EOMMM23A - 0.5*(einsum_kg(H2A.ovoo,r2a,'mcjk,abim->abcijk')...
%                             -einsum_kg(H2A.ovoo,r2a,'mcik,abjm->abcijk')...
%                             -einsum_kg(H2A.ovoo,r2a,'mcji,abkm->abcijk'));
%                         
%      EOMMM23A = EOMMM23A - 0.5*(einsum_kg(I_ovoo,t2a,'mcjk,abim->abcijk')...
%                         -einsum_kg(I_ovoo,t2a,'mcik,abjm->abcijk')...
%                         -einsum_kg(I_ovoo,t2a,'mcji,abkm->abcijk'));
%                     
%      EOMMM23A = EOMMM23A + einsum_kg(I_vvov,t2a,'abie,ecjk->abcijk')...
%                        -einsum_kg(I_vvov,t2a,'abje,ecik->abcijk')...
%                        -einsum_kg(I_vvov,t2a,'abke,ecji->abcijk');
%      
%      % A(abc)
%      EOMMM23A = EOMMM23A + ...
%               -permute(EOMMM23A,[2,1,3,4,5,6])...
%               -permute(EOMMM23A,[1,3,2,4,5,6])...
%               -permute(EOMMM23A,[3,2,1,4,5,6])...
%               +permute(EOMMM23A,[2,3,1,4,5,6])...
%               +permute(EOMMM23A,[3,1,2,4,5,6]);
    
    EOMMM23A = 0.0;
   
    EOMMM23A = -einsum_kg(H3A.vvoooo,r1a,'acmikj,bm->abcijk')...
               +einsum_kg(H3A.vvoooo,r1a,'bcmikj,am->abcijk')... % (ab)
               +einsum_kg(H3A.vvoooo,r1a,'abmikj,cm->abcijk');   % (bc)
                
    EOMMM23A = EOMMM23A +einsum_kg(H3A.vvvoov,r1a,'acbike,ej->abcijk')...
                        -einsum_kg(H3A.vvvoov,r1a,'acbije,ek->abcijk')... % (kj)
                        -einsum_kg(H3A.vvvoov,r1a,'acbjke,ei->abcijk');   % (ij)
              
    EOMMM23A = EOMMM23A +0.5*einsum_kg(H3A.vvvvvo,r2a,'abcefk,efij->abcijk')...
                        -0.5*einsum_kg(H3A.vvvvvo,r2a,'abcefi,efkj->abcijk')... % (ki)
                        -0.5*einsum_kg(H3A.vvvvvo,r2a,'abcefj,efik->abcijk');   % (kj)
            
    EOMMM23A = EOMMM23A +0.5*einsum_kg(H3A.oovooo,r2a,'mncijk,abmn->abcijk')... 
                        -0.5*einsum_kg(H3A.oovooo,r2a,'mnaijk,cbmn->abcijk')... % (ac)
                        -0.5*einsum_kg(H3A.oovooo,r2a,'mnbijk,acmn->abcijk');   % (bc)
            
    EOMMM23A = EOMMM23A +einsum_kg(H3A.vvooov,r2a,'abmije,ecmk->abcijk')... 
                        -einsum_kg(H3A.vvooov,r2a,'cbmije,eamk->abcijk')... % (ac)
                        -einsum_kg(H3A.vvooov,r2a,'acmije,ebmk->abcijk')... % (bc)
                        -einsum_kg(H3A.vvooov,r2a,'abmkje,ecmi->abcijk')... % (ik)
                        -einsum_kg(H3A.vvooov,r2a,'abmike,ecmj->abcijk')... % (jk)
                        +einsum_kg(H3A.vvooov,r2a,'cbmkje,eami->abcijk')... % (ik)(ac)
                        +einsum_kg(H3A.vvooov,r2a,'cbmike,eamj->abcijk')... % (jk)(ac)
                        +einsum_kg(H3A.vvooov,r2a,'acmkje,ebmi->abcijk')... % (ik)(bc)
                        +einsum_kg(H3A.vvooov,r2a,'acmike,ebmj->abcijk');   % (jk)(bc)
                    
     EOMMM23A = EOMMM23A +einsum_kg(H3B.vvooov,r2b,'abmije,cekm->abcijk')... 
                         -einsum_kg(H3B.vvooov,r2b,'cbmije,aekm->abcijk')... % (ac)
                         -einsum_kg(H3B.vvooov,r2b,'acmije,bekm->abcijk')... % (bc)
                         -einsum_kg(H3B.vvooov,r2b,'abmkje,ceim->abcijk')... % (ik)
                         -einsum_kg(H3B.vvooov,r2b,'abmike,cejm->abcijk')... % (jk)
                         +einsum_kg(H3B.vvooov,r2b,'cbmkje,aeim->abcijk')... % (ik)(ac)
                         +einsum_kg(H3B.vvooov,r2b,'cbmike,aejm->abcijk')... % (jk)(ac)
                         +einsum_kg(H3B.vvooov,r2b,'acmkje,beim->abcijk')... % (ik)(bc)
                         +einsum_kg(H3B.vvooov,r2b,'acmike,bejm->abcijk');   % (jk)(bc)
            
     EOMMM23A = EOMMM23A -einsum_kg(H2A.vooo,r2a,'amij,bcmk->abcijk')...
                         +einsum_kg(H2A.vooo,r2a,'amkj,bcmi->abcijk')... % (ki)
                         +einsum_kg(H2A.vooo,r2a,'amik,bcmj->abcijk')... % (kj)
                         +einsum_kg(H2A.vooo,r2a,'bmij,acmk->abcijk')... % (ab)
                         +einsum_kg(H2A.vooo,r2a,'cmij,bamk->abcijk')... % (ac)
                         -einsum_kg(H2A.vooo,r2a,'bmkj,acmi->abcijk')... % (ki)(ab)
                         -einsum_kg(H2A.vooo,r2a,'cmkj,bami->abcijk')... % (ki)(ac)
                         -einsum_kg(H2A.vooo,r2a,'bmik,acmj->abcijk')... % (kj)(ab)
                         -einsum_kg(H2A.vooo,r2a,'cmik,bamj->abcijk');   % (ij)(ac)
             
     EOMMM23A = EOMMM23A +einsum_kg(H2A.vvov,r2a,'abie,ecjk->abcijk')...
                         -einsum_kg(H2A.vvov,r2a,'cbie,eajk->abcijk')... % (ac)
                         -einsum_kg(H2A.vvov,r2a,'acie,ebjk->abcijk')... % (bc)
                         -einsum_kg(H2A.vvov,r2a,'abke,ecji->abcijk')... % (ik)
                         -einsum_kg(H2A.vvov,r2a,'abje,ecik->abcijk')... % (ij)
                         +einsum_kg(H2A.vvov,r2a,'cbke,eaji->abcijk')... % (ac)(ik)
                         +einsum_kg(H2A.vvov,r2a,'cbje,eaik->abcijk')... % (ac)(ij)
                         +einsum_kg(H2A.vvov,r2a,'acke,ebji->abcijk')... % (bc)(ik)
                         +einsum_kg(H2A.vvov,r2a,'acje,ebik->abcijk');   % (bc)(ij)

     I1 = einsum_kg(VA,r1a,'mnef,fn->me');
     EOMMM23A = EOMMM23A - einsum_kg(einsum_kg(I1,t2a,'me,ecjk->mcjk'),t2a,'mcjk,abim->abcijk') ...
                         + einsum_kg(einsum_kg(I1,t2a,'me,ecik->mcik'),t2a,'mcik,abjm->abcijk') ... % (ij) 
                         + einsum_kg(einsum_kg(I1,t2a,'me,ecji->mcji'),t2a,'mcji,abkm->abcijk') ... % (ik), only need A(i/jk) i think...
                         + einsum_kg(einsum_kg(I1,t2a,'me,eajk->majk'),t2a,'majk,cbim->abcijk') ... % (ac)
                         + einsum_kg(einsum_kg(I1,t2a,'me,ebjk->mbjk'),t2a,'mbjk,acim->abcijk') ... % (bc)
                         - einsum_kg(einsum_kg(I1,t2a,'me,eaik->maik'),t2a,'maik,cbjm->abcijk') ... % (ij)(ac)
                         - einsum_kg(einsum_kg(I1,t2a,'me,ebik->mbik'),t2a,'mbik,acjm->abcijk') ... % (ij)(bc)
                         - einsum_kg(einsum_kg(I1,t2a,'me,eaji->maji'),t2a,'maji,cbkm->abcijk') ... % (ik)(ac)
                         - einsum_kg(einsum_kg(I1,t2a,'me,ebji->mbji'),t2a,'mbji,ackm->abcijk');    % (ik)(bc)
     
     I2 = einsum_kg(VB,r1b,'mnef,fn->me');
     EOMMM23A = EOMMM23A - einsum_kg(einsum_kg(I2,t2a,'me,ecjk->mcjk'),t2a,'mcjk,abim->abcijk') ...
                         + einsum_kg(einsum_kg(I2,t2a,'me,ecik->mcik'),t2a,'mcik,abjm->abcijk') ... % (ij) 
                         + einsum_kg(einsum_kg(I2,t2a,'me,ecji->mcji'),t2a,'mcji,abkm->abcijk') ... % (ik), only need A(i/jk) i think...
                         + einsum_kg(einsum_kg(I2,t2a,'me,eajk->majk'),t2a,'majk,cbim->abcijk') ... % (ac)
                         + einsum_kg(einsum_kg(I2,t2a,'me,ebjk->mbjk'),t2a,'mbjk,acim->abcijk') ... % (bc)
                         - einsum_kg(einsum_kg(I2,t2a,'me,eaik->maik'),t2a,'maik,cbjm->abcijk') ... % (ij)(ac)
                         - einsum_kg(einsum_kg(I2,t2a,'me,ebik->mbik'),t2a,'mbik,acjm->abcijk') ... % (ij)(bc)
                         - einsum_kg(einsum_kg(I2,t2a,'me,eaji->maji'),t2a,'maji,cbkm->abcijk') ... % (ik)(ac)
                         - einsum_kg(einsum_kg(I2,t2a,'me,ebji->mbji'),t2a,'mbji,ackm->abcijk');    % (ik)(bc)
     
     %EOMMM23 = EOMMM23 + r0*MM23
     
     fprintf(' finished in %4.2f s\n',toc)
             
