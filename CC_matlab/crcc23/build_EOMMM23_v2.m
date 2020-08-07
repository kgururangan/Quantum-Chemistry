function [EOMMM23] = build_EOMMM23_v2(t1,t2,r1,r2,HBar,sys)

%%% DIFFERENCE BETWEEN build_EOMMM23_v2 AND build_EOMMM23 IS THAT THIS
%%% VERSION (COPIED FROM GOUR, WLOCH, PIECUCH PAPER SEEMS TO INCLUDE THE
%%% ANTISYMMETRIZERS A(i/jk)A(c/ab) ON THE 4-BODY HBAR TERM CONTRACTED WITH
%%% R1.

%     H3 = HBar{3}; H2 = HBar{2}; H1 = HBar{1};
% %     H_vvoooo = H3{2,2,1,1,1,1};
% %     H_vvooov = H3{2,2,1,1,1,2};
% %     H_vvvovo = H3{2,2,2,1,2,1};
% %     H_vvvvvo = H3{2,2,2,2,2,1};
% %     H_oovooo = H3{1,1,2,1,1,1};
%     H_oooo = H2{1,1,1,1};
%     H_vvvv = H2{2,2,2,2};
%     H_vvov = H2{2,2,1,2};
%     H_ooov = H2{1,1,1,2};
%     H_vovv = H2{2,1,2,2};
%     H_ovov = H2{1,2,1,2};
%     H_ovoo = H2{1,2,1,1};

    H2 = HBar{2};
    
    
    fprintf('EOM-MMCC(2,3) construction... ')
    
    tic
    
%      H_ov = sys.fov + einsum_kg(sys.Voovv,t1,'imae,em->ia');
    

     
%      V_ovov = sys.Vovov + einsum_kg(sys.Vovvv,t1,'iaeb,ej->iajb');
%      
%      H_vovv = sys.Vvovv - einsum_kg(sys.Voovv,t1,'mibc,am->aibc');
%      
%      H_ooov = sys.Vooov + einsum_kg(sys.Voovv,t1,'ijea,ek->ijka');
     
%      H_vvvv = sys.Vvvvv + 0.5*einsum_kg(sys.Voovv,t2,'mncd,abmn->abcd')...
%               -einsum_kg(H_vovv,t1,'amcd,bm->abcd') + einsum_kg(sys.Vvovv,t1,'bmcd,am->abcd');
%           
%      H_oooo = sys.Voooo + 0.5*einsum_kg(sys.Voovv,t2,'ijef,efkl->ijkl')...
%               -einsum_kg(H_ooov,t1,'ijle,ek->ijkl') + einsum_kg(sys.Vooov,t1,'ijke,el->ijkl');
%           
%      H_ovov = V_ovov - einsum_kg(sys.Voovv,t2,'imeb,eajm->iajb') - einsum_kg(H_ooov,t1,'imjb,am->iajb');
     
%      H_vvov = sys.Vvvov + einsum_kg(sys.Vvvvv,t1,'abec,ei->abic') - einsum_kg(H_ovov,t1,'mbic,am->abic')...
%               +einsum_kg(V_ovov,t1,'maic,bm->abic') - einsum_kg(H_ov,t2,'mc,abim->abic')...
%               +einsum_kg(H_vovv,t2,'bmce,aeim->abic') - einsum_kg(sys.Vvovv,t2,'amce,beim->abic')...
%               +0.5*einsum_kg(H_ooov,t2,'nmic,abnm->abic');
%           
%      H_ovoo = sys.Vovoo + einsum_kg(H_oooo,t1,'mijk,am->iajk') - einsum_kg(sys.Vovov,t1,'iake,ej->iajk')...
%               +einsum_kg(H_ooov,t2,'imje,aekm->iajk') - einsum_kg(H_ooov,t2,'imke,aejm->iajk')...
%               +einsum_kg(H_ov,t2,'ie,eajk->iajk') + einsum_kg(V_ovov,t1,'iaje,ek->iajk')...
%               -0.5*einsum_kg(sys.Vvovv,t2,'aief,efjk->iajk');

     I_ov = -einsum_kg(sys.Voovv,r1,'miae,em->ia');
     
     I_ovoo = einsum_kg(I_ov,t2,'ie,aekj->iajk')...
              +0.5*einsum_kg(H2{2,1,2,2},r2,'aief,efkj->iajk')...
              -einsum_kg(H2{1,1,1,1},r1,'imjk,am->iajk')...
              +einsum_kg(H2{1,1,1,2},r2,'imje,aekm->iajk')+einsum_kg(H2{1,2,1,2},r1,'iaje,ek->iajk')...
              -einsum_kg(H2{1,1,1,2},r2,'imke,aejm->iajk')-einsum_kg(H2{1,2,1,2},r1,'iake,ej->iajk');
     
     I_vvov = 0.5*einsum_kg(H2{2,2,2,2},r1,'abec,ei->abic')...
              -einsum_kg(H2{2,1,2,2},r2,'bmec,aeim->abic')...
              +0.25*einsum_kg(H2{1,1,1,2},r2,'mnic,abmn->abic')...
              -einsum_kg(H2{1,2,1,2},r1,'mbic,am->abic');
          
     EOMMM23 = 0.5*(einsum_kg(H2{2,2,1,2},r2,'abie,ecjk->abcijk')...
                    -einsum_kg(H2{2,2,1,2},r2,'abje,ecik->abcijk')...
                    -einsum_kg(H2{2,2,1,2},r2,'abke,ecji->abcijk'));
     
     EOMMM23 = EOMMM23 - 0.5*(einsum_kg(H2{1,2,1,1},r2,'mcjk,abim->abcijk')...
                            -einsum_kg(H2{1,2,1,1},r2,'mcik,abjm->abcijk')...
                            -einsum_kg(H2{1,2,1,1},r2,'mcji,abkm->abcijk'));
                        
     EOMMM23 = EOMMM23 - 0.5*(einsum_kg(I_ovoo,t2,'mcjk,abim->abcijk')...
                        -einsum_kg(I_ovoo,t2,'mcik,abjm->abcijk')...
                        -einsum_kg(I_ovoo,t2,'mcji,abkm->abcijk'));
                    
     EOMMM23 = EOMMM23 + einsum_kg(I_vvov,t2,'abie,ecjk->abcijk')...
                       -einsum_kg(I_vvov,t2,'abje,ecik->abcijk')...
                       -einsum_kg(I_vvov,t2,'abke,ecji->abcijk');
     
     % A(abc)
     EOMMM23 = EOMMM23 + ...
              -permute(EOMMM23,[2,1,3,4,5,6])...
              -permute(EOMMM23,[1,3,2,4,5,6])...
              -permute(EOMMM23,[3,2,1,4,5,6])...
              +permute(EOMMM23,[2,3,1,4,5,6])...
              +permute(EOMMM23,[3,1,2,4,5,6]);
     
     fprintf(' finished in %4.2f s\n',toc)
             

end

