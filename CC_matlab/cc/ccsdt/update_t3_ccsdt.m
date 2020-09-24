function [t3] = update_t3_ccsdt(t1,t2,t3,sys)
    
    % Hbar diagram intermediates
    h_me = sys.fov + einsum_kg(sys.Voovv,t1,'mnef,fn->me');
    
    tau = t2 + 0.5*(einsum_kg(t1,t1,'ai,bj->abij') - ...
                    einsum_kg(t1,t1,'aj,bi->abij') - ...
                    einsum_kg(t1,t1,'bi,aj->abij') + ...
                    einsum_kg(t1,t1,'bj,ai->abij')); 
                
    Wbar_mnij = sys.Voooo + ...
                einsum_kg(sys.Vooov,t1,'mnie,ej->mnij') - einsum_kg(sys.Vooov,t1,'mnje,ei->mnij') + ...
                0.25*einsum_kg(sys.Voovv,tau,'mnef,efij->mnij');
            
    Wbar_abef = sys.Vvvvv - ...
                einsum_kg(sys.Vvovv,t1,'amef,bm->abef') + einsum_kg(sys.Vvovv,t1,'bmef,am->abef') + ...
                0.25*einsum_kg(sys.Voovv,tau,'mnef,abmn->abef');
    
    h_mnij = Wbar_mnij + 0.25*einsum_kg(sys.Voovv,tau,'mnef,efij->mnij'); 
    
    h_abef = Wbar_abef + 0.25*einsum_kg(sys.Voovv,tau,'mnef,abmn->abef'); 
    
    h_mbij =   sys.Vovoo - einsum_kg(h_me,t2,'me,beij->mbij') - einsum_kg(h_mnij,t1,'mnij,bn->mbij') + ...
                    0.5*einsum_kg(sys.Vovvv,tau,'mbef,efij->mbij') + ...
                    einsum_kg(sys.Vooov,t2,'mnie,bejn->mbij') - einsum_kg(sys.Vooov,t2,'mnje,bein->mbij') + ...
                    einsum_kg(sys.Vovvo,t1,'mbej,ei->mbij') - einsum_kg(sys.Vovvo,t1,'mbei,ej->mbij') - ...
                    einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,bfnj->mbej'),t1,'mbej,ei->mbij') + ...
                    einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,bfni->mbei'),t1,'mbei,ej->mbij'); % 
                
    h_amij = -permute(h_mbij,[2,1,3,4]);
                
    h_abej =   sys.Vvvvo - einsum_kg(h_me,t2,'me,abmi->abei') + einsum_kg(h_abef,t1,'abef,fi->abei') + ...
                    0.5*einsum_kg(sys.Voovo,tau,'mnei,abmn->abei') - ...
                    einsum_kg(sys.Vovvv,t2,'mbef,afmi->abei') + einsum_kg(sys.Vovvv,t2,'maef,bfmi->abei') - ...
                    einsum_kg(sys.Vovvo,t1,'mbei,am->abei') + einsum_kg(sys.Vovvo,t1,'maei,bm->abei') + ...
                    einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,bfni->mbei'),t1,'mbei,am->abei') - ...
                    einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,afni->maei'),t1,'maei,bm->abei'); 
                
    h_abie = -permute(h_abej,[1,2,4,3]); 
        
    % intermediates for projection onto T3
    I1 = sys.Vvoov ...
         - einsum_kg(sys.Voovo,t1,'mnei,an->amie') ...
         + einsum_kg(sys.Vvovv,t1,'amfe,fi->amie') ...
         + einsum_kg(sys.Voovv,t2,'mnef,afin->amie') ...
         - einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fi->mnei'),t1,'mnei,an->amie');
    
    I2 = sys.Voooo ...
         + einsum_kg(sys.Vooov,t1,'mnie,ej->mnij') - einsum_kg(sys.Vooov,t1,'mnje,ei->mnij') ...
         + 0.5*einsum_kg(sys.Voovv,t2,'mnef,efij->mnij') ...
         + einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,ei->mnif'),t1,'mnif,fj->mnij');
     
    I3 = sys.Vvvvv ...
         + einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,am->anef'),t1,'anef,bn->abef') ...
         + 0.5*einsum_kg(sys.Voovv,t2,'mnef,abmn->abef') ...
         - einsum_kg(sys.Vvovv,t1,'amef,bm->abef') + einsum_kg(sys.Vvovv,t1,'bmef,am->abef');
     
    I4 = einsum_kg(sys.Voovv,t2,'mnfe,afij->amnije'); % W_amnije
    
    I5 = -einsum_kg(sys.Voovv,t2,'mnfe,abin->abmief'); % W_abmief
    
    I6 = sys.foo_masked ...
         + einsum_kg(sys.fov,t1,'me,ei->mi') ...
         + einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fn->me'),t1,'me,ei->mi') ...
         + 0.5*einsum_kg(sys.Voovv,t2,'mnef,efin->mi') ...
         + einsum_kg(sys.Vooov,t1,'mnif,fn->mi');
    
    I7 = sys.fvv_masked ...
         - einsum_kg(sys.fov,t1,'me,am->ae') ...
         - einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fn->me'),t1,'me,am->ae') ...
         - 0.5*einsum_kg(sys.Voovv,t2,'mnef,afmn->ae') ...
         + einsum_kg(sys.Vvovv,t1,'amef,fm->ae');
     
    I12 = h_abie-einsum_kg(h_me,t2,'me,abim->abie');
    
    % build T3
    
    % < phi_{ijkabc} | H_{CCSD} | 0 >
    % = -A(k/ij)A(a/bc) h(amij)*t(bcmk) + A(i/jk)A(c/ab)(h(abie)-h(me)*t(abim))*t(ecjk)
%     MM23 =   -einsum_kg(h_amij,t2,'amij,bcmk->abcijk') ... % (1)
%              +einsum_kg(h_amij,t2,'amkj,bcmi->abcijk') ... % (ik)
%              +einsum_kg(h_amij,t2,'amik,bcmj->abcijk') ... % (jk)
%              +einsum_kg(h_amij,t2,'cmij,bamk->abcijk') ... % (ac)
%              +einsum_kg(h_amij,t2,'bmij,acmk->abcijk') ... % (ab)
%              -einsum_kg(h_amij,t2,'bmkj,acmi->abcijk') ... % (ab)(ik)
%              -einsum_kg(h_amij,t2,'cmkj,bami->abcijk') ... % (ac)(ik)
%              -einsum_kg(h_amij,t2,'bmik,acmj->abcijk') ... % (ab)(jk)
%              -einsum_kg(h_amij,t2,'cmik,bamj->abcijk');    % (ac)(jk)
%          
%     MM23 =  MM23 +einsum_kg(I12,t2,'abie,ecjk->abcijk') ... % (1)
%                  -einsum_kg(I12,t2,'abje,ecik->abcijk') ... % (ij)
%                  -einsum_kg(I12,t2,'abke,ecji->abcijk') ... % (ik)
%                  -einsum_kg(I12,t2,'cbie,eajk->abcijk') ... % (ac)
%                  -einsum_kg(I12,t2,'acie,ebjk->abcijk') ... % (bc)
%                  +einsum_kg(I12,t2,'cbje,eaik->abcijk') ... % (ac)(ij)
%                  +einsum_kg(I12,t2,'acje,ebik->abcijk') ... % (bc)(ij)
%                  +einsum_kg(I12,t2,'cbke,eaji->abcijk') ... % (ac)(ik)
%                  +einsum_kg(I12,t2,'acke,ebji->abcijk');    % (bc)(ik)
    
%     D1 = -einsum_kg(sys.Vvooo,t2,'amij,bcmk->abcijk');
%     
%     D2 = einsum_kg(sys.Vvvvo,t2,'bcek,aeij->abcijk');
%     
%     D3 = -einsum_kg(einsum_kg(sys.Vvoov,t1,'amie,ej->amij'),t2,'amij,bcmk->abcijk');
%     D3 = D3 - permute(D3,[1,2,3,5,4,6]);
%     
%     D4 = -einsum_kg(einsum_kg(sys.Vovvo,t1,'mcek,bm->bcek'),t2,'bcek,aeij->abcijk');
%     D4 = D4 - permute(D4,[1,3,2,4,5,6]);
%     
%     D5 = einsum_kg(einsum_kg(sys.Vvvvv,t1,'bcef,fk->bcek'),t2,'bcek,aeij->abcijk');
%     
%     D6 = einsum_kg(einsum_kg(sys.Voooo,t1,'nmij,an->amij'),t2,'amij,bcmk->abcijk');
%     
%     D7 = -einsum_kg(einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fi->amie'),t1,'amie,ej->amij'),t2,'amij,bcmk->abcijk');
%     
%     D8 = einsum_kg(einsum_kg(einsum_kg(sys.Vooov,t1,'nmie,ej->nmij'),t1,'nmij,an->amij'),t2,'amij,bcmk->abcijk');
%     
%     D9 = einsum_kg(einsum_kg(einsum_kg(sys.Voovo,t1,'mnek,cn->mcek'),t1,'mcek,bm->bcek'),t2,'bcek,aeij->abcijk');
%     
%     D10 = -einsum_kg(einsum_kg(einsum_kg(sys.Vovvv,t1,'mcef,fk->mcek'),t1,'mcek,bm->bcek'),t2,'bcek,aeij->abcijk');
%     
%     D11 = -einsum_kg(einsum_kg(einsum_kg(sys.Vvovv,t1,'bnef,cn->bcef'),t1,'bcef,fk->bcek'),t2,'bcek,aeij->abcijk');
%     
%     D12 = einsum_kg(einsum_kg(einsum_kg(sys.Voovo,t1,'nmfj,an->amfj'),t1,'amfj,fi->amij'),t2,'amij,bcmk->abcijk');
%     
%     D13 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,cn->mcef'),t1,'mcef,bm->bcef'),t1,'bcef,fk->bcek'),t2,'bcek,aeij->abcijk');
%     
%     D14 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,ej->nmfj'),t1,'nmfj,an->amfj'),t1,'amfj,fi->amij'),t2,'amij,bcmk->abcijk');
%     
%     D15 = -einsum_kg(einsum_kg(sys.fov,t2,'me,bcmk->bcek'),t2,'bcek,aeij->abcijk');
%     
%     D16 = 0.5*einsum_kg(einsum_kg(sys.Voovo,t2,'mnek,bcmn->bcek'),t2,'bcek,aeij->abcijk');
%     
%     D17 = -0.5*einsum_kg(einsum_kg(sys.Vvovv,t2,'amfe,feij->amij'),t2,'amij,bcmk->abcijk');
%     
%     D18 = -einsum_kg(einsum_kg(sys.Vooov,t2,'mnjf,afin->amij'),t2,'amij,bcmk->abcijk');
%     D18 = D18 - permute(D18,[1,2,3,5,4,6]);
%     
%     D19 = einsum_kg(einsum_kg(sys.Vvovv,t2,'bnef,fcnk->bcek'),t2,'bcek,aeij->abcijk');
%     D19 = D19 - permute(D19,[1,3,2,4,5,6]);
%     
%     D20 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,bm->bnef'),t1,'bnef,cn->bcef'),t1,'bcef,fk->bcek'),t2,'bcek,aeij->abcijk');
%     
%     D21 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,fi->nmie'),t1,'nmie,an->amie'),t1,'amie,ej->amij'),t2,'amij,bcmk->abcijk');
%     
%     D22 = -einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fn->me'),t2,'me,bcmk->bcek'),t2,'bcek,aeij->abcijk');
% 
%     % A(k/ij)A(a/bc)
%     MM23 = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14 + D15 + D16 + D17 + D18 + D19 + D20 + D21 + D22;
%     MM23 = MM23 - permute(MM23,[2,1,3,4,5,6]) - permute(MM23,[3,2,1,4,5,6]);
%     MM23 = MM23 - permute(MM23,[1,2,3,6,5,4]) - permute(MM23,[1,2,3,4,6,5]);

    D1 = sys.Vvooo;
    %D1 = -einsum_kg(sys.Vvooo,t2,'amij,bcmk->abcijk');
    
    D2 = sys.Vvvvo;
    %D2 = einsum_kg(sys.Vvvvo,t2,'bcek,aeij->abcijk');
    
    D3 = einsum_kg(sys.Vvoov,t1,'amie,ej->amij') - einsum_kg(sys.Vvoov,t1,'amje,ei->amij');
    %D3 = -einsum_kg(D3,t2,'amij,bcmk->abcijk');
    
    D4 = -einsum_kg(sys.Vovvo,t1,'mcek,bm->bcek') + einsum_kg(sys.Vovvo,t1,'mbek,cm->bcek');
    %D4 = einsum_kg(D4,t2,'bcek,aeij->abcijk');
    
    D5 = einsum_kg(sys.Vvvvv,t1,'bcef,fk->bcek');
    %D5 = einsum_kg(einsum_kg(sys.Vvvvv,t1,'bcef,fk->bcek'),...
    %                                   t2,'bcek,aeij->abcijk');
    
    D6 = -einsum_kg(sys.Voooo,t1,'nmij,an->amij');
    %D6 = einsum_kg(einsum_kg(sys.Voooo,t1,'nmij,an->amij'),...
    %                                   t2,'amij,bcmk->abcijk');
    
    D7 = einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fi->amie'),t1,'amie,ej->amij');
    %D7 = -einsum_kg(einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fi->amie'),...
    %                                              t1,'amie,ej->amij'),...
    %                                              t2,'amij,bcmk->abcijk');
    
    D8 = -einsum_kg(einsum_kg(sys.Vooov,t1,'nmie,ej->nmij'),t1,'nmij,an->amij') + einsum_kg(einsum_kg(sys.Vooov,t1,'nmje,ei->nmij'),t1,'nmij,an->amij');
    %D8 = -einsum_kg(D8,t2,'amij,bcmk->abcijk');                                     
    
    D9 = einsum_kg(einsum_kg(sys.Voovo,t1,'mnek,cn->mcek'),t1,'mcek,bm->bcek');
    %D9 = einsum_kg(einsum_kg(einsum_kg(sys.Voovo,t1,'mnek,cn->mcek'),...
    %                                             t1,'mcek,bm->bcek'),...
    %                                             t2,'bcek,aeij->abcijk');
    
    D10 = -einsum_kg(einsum_kg(sys.Vovvv,t1,'mcef,fk->mcek'),t1,'mcek,bm->bcek') + einsum_kg(einsum_kg(sys.Vovvv,t1,'mbef,fk->mbek'),t1,'mbek,cm->bcek');                                    
    %D10 = einsum_kg(D10,t2,'bcek,aeij->abcijk');                                       
           
    D13 = einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,cn->mcef'),t1,'mcef,bm->bcef'),t1,'bcef,fk->bcek');
    %D13 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,cn->mcef'),...
    %                                                        t1,'mcef,bm->bcef'),...
    %                                                        t1,'bcef,fk->bcek'),...
    %                                                        t2,'bcek,aeij->abcijk');
    
    D14 = -einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,ej->nmfj'),...
                                                   t1,'nmfj,an->amfj'),...
                                                   t1,'amfj,fi->amij');
    %D14 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,ej->nmfj'),...
    %                                                        t1,'nmfj,an->amfj'),...
    %                                                        t1,'amfj,fi->amij'),...
    %                                                        t2,'amij,bcmk->abcijk');
    
    D15 = -einsum_kg(sys.fov,t2,'me,bcmk->bcek');
    %D15 = -einsum_kg(einsum_kg(sys.fov,t2,'me,bcmk->bcek'),...
    %                                   t2,'bcek,aeij->abcijk');
    
    D16 = 0.5*einsum_kg(sys.Voovo,t2,'mnek,bcmn->bcek');
    %D16 = 0.5*einsum_kg(einsum_kg(sys.Voovo,t2,'mnek,bcmn->bcek'),...
    %                                        t2,'bcek,aeij->abcijk');
    
    D17 = 0.5*einsum_kg(sys.Vvovv,t2,'amfe,feij->amij');
    %D17 = -0.5*einsum_kg(einsum_kg(sys.Vvovv,t2,'amfe,feij->amij'),...
    %                                         t2,'amij,bcmk->abcijk');
    
    D18 = einsum_kg(sys.Vooov,t2,'mnjf,afin->amij') - einsum_kg(sys.Vooov,t2,'mnif,afjn->amij');                    
    %D18 = -einsum_kg(D18,t2,'amij,bcmk->abcijk');
    
    D19 = einsum_kg(sys.Vvovv,t2,'bnef,fcnk->bcek') - einsum_kg(sys.Vvovv,t2,'cnef,fbnk->bcek');
    %D19 = einsum_kg(D19,t2,'bcek,aeij->abcijk');
    
    D22 = -einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fn->me'),...
                                         t2,'me,bcmk->bcek');
    %D22 = -einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fn->me'),...
    %                                               t2,'me,bcmk->bcek'),...
    %                                               t2,'bcek,aeij->abcijk');
                                               
                                              
    D23 = -0.5*einsum_kg(einsum_kg(sys.Voovv,t2,'nmfe,efji->nmij'),t1,'nmij,an->amij');    
    %D23 = -einsum_kg(D23,t2,'amij,bcmk->abcijk');
    
    D24 = 0.5*einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,bcmn->bcef'),t1,'bcef,fk->bcek');
    %D24 = einsum_kg(D24,t2,'bcek,aeij->abcijk');
    
    D25 = einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,afin->amie'),t1,'amie,ej->amij') ...
                -einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,afjn->amje'),t1,'amje,ei->amij');
    %D25 = -einsum_kg(D25,t2,'amij,bcmk->abcijk');
            
    D26 = -einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,fcnk->mcek'),t1,'mcek,bm->bcek') ...
                + einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,fbnk->mbek'),t1,'mbek,cm->bcek');
    %D26 = einsum_kg(D26,t2,'bcek,aeij->abcijk');
    
    h_vooo = D1 + D3 + D6 + D7 +  D8 + D14 + D17 + D18 + D23 + D25;
    h_vvvo = D2 + D4 + D5 + D9 + D10 + D13 + D16 + D19 + D24 + D26; 
    %h_ov_t2 = D15 + D22;
    h_ov = sys.fov + einsum_kg(sys.Voovv,t1,'mnef,fn->me');

    h_vooo_t2 = -einsum_kg(h_vooo,t2,'amij,bcmk->abcijk');
    h_vvvo_t2 = einsum_kg(h_vvvo,t2,'bcek,aeij->abcijk');
    h_ov_t2_t2 = -einsum_kg(h_ov,t2,'me,bcmk->bcek');
    h_ov_t2_t2 = einsum_kg(h_ov_t2_t2,t2,'bcek,aeij->abcijk');
    
    h_vooo_t2 = h_vooo_t2 - permute(h_vooo_t2,[2,1,3,4,5,6]) - permute(h_vooo_t2,[3,2,1,4,5,6]);
    h_vooo_t2 = h_vooo_t2 - permute(h_vooo_t2,[1,2,3,6,5,4]) - permute(h_vooo_t2,[1,2,3,4,6,5]);
    
    h_vvvo_t2 = h_vvvo_t2 - permute(h_vvvo_t2,[2,1,3,4,5,6]) - permute(h_vvvo_t2,[3,2,1,4,5,6]);
    h_vvvo_t2 = h_vvvo_t2 - permute(h_vvvo_t2,[1,2,3,6,5,4]) - permute(h_vvvo_t2,[1,2,3,4,6,5]);
    
    h_ov_t2_t2 = h_ov_t2_t2 - permute(h_ov_t2_t2,[2,1,3,4,5,6]) - permute(h_ov_t2_t2,[3,2,1,4,5,6]);
    h_ov_t2_t2 = h_ov_t2_t2 - permute(h_ov_t2_t2,[1,2,3,6,5,4]) - permute(h_ov_t2_t2,[1,2,3,4,6,5]);

    MM23 = h_vooo_t2 + h_vvvo_t2 + h_ov_t2_t2;
    
    % < phi_{ijkabc} | [H_N(T3 + T1*T3 + T1^2*T3 + T2*T3)]_C | 0 >
    
    % part 1
    % A(c/ab)[ A(k/ij)I1(cmke)t3(abeijm) + 0.5*I3(abef)t3(efcijk) 
    %          + 0.5*A(i/jk)*I5(abmief)t3(fecmjk) + I7(ce)t3(abeijk) ]
    TEMP3_abc = einsum_kg(I1,t3,'cmke,abeijm->abcijk') - einsum_kg(I1,t3,'cmie,abekjm->abcijk') - einsum_kg(I1,t3,'cmje,abeikm->abcijk') ...
                + 0.5*einsum_kg(I3,t3,'abef,efcijk->abcijk') ...
                + 0.5*einsum_kg(I5,t3,'abmief,fecmjk->abcijk') - 0.5*einsum_kg(I5,t3,'abmjef,fecmik->abcijk') - 0.5*einsum_kg(I5,t3,'abmkef,fecmji->abcijk') ...
                + einsum_kg(I7,t3,'ce,abeijk->abcijk');
            
    TEMP3_cba = einsum_kg(I1,t3,'amke,cbeijm->abcijk') - einsum_kg(I1,t3,'amie,cbekjm->abcijk') - einsum_kg(I1,t3,'amje,cbeikm->abcijk') ...
                + 0.5*einsum_kg(I3,t3,'cbef,efaijk->abcijk') ...
                + 0.5*einsum_kg(I5,t3,'cbmief,feamjk->abcijk') - 0.5*einsum_kg(I5,t3,'cbmjef,feamik->abcijk') - 0.5*einsum_kg(I5,t3,'cbmkef,feamji->abcijk') ...
                + einsum_kg(I7,t3,'ae,cbeijk->abcijk');
    
    TEMP3_acb = einsum_kg(I1,t3,'bmke,aceijm->abcijk') - einsum_kg(I1,t3,'bmie,acekjm->abcijk') - einsum_kg(I1,t3,'bmje,aceikm->abcijk') ...
                + 0.5*einsum_kg(I3,t3,'acef,efbijk->abcijk') ...
                + 0.5*einsum_kg(I5,t3,'acmief,febmjk->abcijk') - 0.5*einsum_kg(I5,t3,'acmjef,febmik->abcijk') - 0.5*einsum_kg(I5,t3,'acmkef,febmji->abcijk') ...
                + einsum_kg(I7,t3,'be,aceijk->abcijk');
            
    % part 2
    % A(k/ij)[ 0.5*I2(mnij)t3(abcmnk) - 0.5*A(a/bc)*I4(amnije)t3(ebcnmk)
    %          - I6(mk)t3(abcijm) ]
    TEMP3_ijk = 0.5*einsum_kg(I2,t3,'mnij,abcmnk->abcijk') ...
                - 0.5*einsum_kg(I4,t3,'amnije,ebcnmk->abcijk') + 0.5*einsum_kg(I4,t3,'bmnije,eacnmk->abcijk') + 0.5*einsum_kg(I4,t3,'cmnije,ebanmk->abcijk') ... 
                - einsum_kg(I6,t3,'mk,abcijm->abcijk');
            
    TEMP3_kji = 0.5*einsum_kg(I2,t3,'mnkj,abcmni->abcijk') ...
                - 0.5*einsum_kg(I4,t3,'amnkje,ebcnmi->abcijk') + 0.5*einsum_kg(I4,t3,'bmnkje,eacnmi->abcijk') + 0.5*einsum_kg(I4,t3,'cmnkje,ebanmi->abcijk') ... 
                - einsum_kg(I6,t3,'mi,abckjm->abcijk');
            
    TEMP3_ikj = 0.5*einsum_kg(I2,t3,'mnik,abcmnj->abcijk') ...
                - 0.5*einsum_kg(I4,t3,'amnike,ebcnmj->abcijk') + 0.5*einsum_kg(I4,t3,'bmnike,eacnmj->abcijk') + 0.5*einsum_kg(I4,t3,'cmnike,ebanmj->abcijk') ... 
                - einsum_kg(I6,t3,'mj,abcikm->abcijk');
            
    TEMP3_3 = TEMP3_abc - TEMP3_cba - TEMP3_acb + TEMP3_ijk - TEMP3_kji - TEMP3_ikj;
            
    X_abcijk = MM23 + TEMP3_3;
    
    for a = 1:sys.Nunocc
        for b = a+1:sys.Nunocc
            for c = b+1:sys.Nunocc
                for i = 1:sys.Nocc
                    for j = i+1:sys.Nocc
                        for k = j+1:sys.Nocc
                            
                            % (1)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(a,b,c,i,j,k) = X_abcijk(a,b,c,i,j,k)/...
                                (sys.foo(i,i)+sys.foo(j,j)+sys.foo(k,k)-sys.fvv(a,a)-sys.fvv(b,b)-sys.fvv(c,c));
                            t3(a,b,c,k,i,j) = t3(a,b,c,i,j,k);
                            t3(a,b,c,j,k,i) = t3(a,b,c,i,j,k);
                            t3(a,b,c,i,k,j) = -t3(a,b,c,i,j,k);
                            t3(a,b,c,j,i,k) = -t3(a,b,c,i,j,k);
                            t3(a,b,c,k,j,i) = -t3(a,b,c,i,j,k);
                            
                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(b,a,c,i,j,k) = -t3(a,b,c,i,j,k);
                            t3(b,a,c,k,i,j) = -t3(a,b,c,i,j,k);
                            t3(b,a,c,j,k,i) = -t3(a,b,c,i,j,k);
                            t3(b,a,c,i,k,j) = t3(a,b,c,i,j,k);
                            t3(b,a,c,j,i,k) = t3(a,b,c,i,j,k);
                            t3(b,a,c,k,j,i) = t3(a,b,c,i,j,k);
                            
                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(a,c,b,i,j,k) = -t3(a,b,c,i,j,k);
                            t3(a,c,b,k,i,j) = -t3(a,b,c,i,j,k);
                            t3(a,c,b,j,k,i) = -t3(a,b,c,i,j,k);
                            t3(a,c,b,i,k,j) = t3(a,b,c,i,j,k);
                            t3(a,c,b,j,i,k) = t3(a,b,c,i,j,k);
                            t3(a,c,b,k,j,i) = t3(a,b,c,i,j,k);
                            
                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(c,b,a,i,j,k) = -t3(a,b,c,i,j,k);
                            t3(c,b,a,k,i,j) = -t3(a,b,c,i,j,k);
                            t3(c,b,a,j,k,i) = -t3(a,b,c,i,j,k);
                            t3(c,b,a,i,k,j) = t3(a,b,c,i,j,k);
                            t3(c,b,a,j,i,k) = t3(a,b,c,i,j,k);
                            t3(c,b,a,k,j,i) = t3(a,b,c,i,j,k);
                            
                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(b,c,a,i,j,k) = t3(a,b,c,i,j,k);
                            t3(b,c,a,k,i,j) = t3(a,b,c,i,j,k);
                            t3(b,c,a,j,k,i) = t3(a,b,c,i,j,k);
                            t3(b,c,a,i,k,j) = -t3(a,b,c,i,j,k);
                            t3(b,c,a,j,i,k) = -t3(a,b,c,i,j,k);
                            t3(b,c,a,k,j,i) = -t3(a,b,c,i,j,k);
                            
                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(c,a,b,i,j,k) = t3(a,b,c,i,j,k);
                            t3(c,a,b,k,i,j) = t3(a,b,c,i,j,k);
                            t3(c,a,b,j,k,i) = t3(a,b,c,i,j,k);
                            t3(c,a,b,i,k,j) = -t3(a,b,c,i,j,k);
                            t3(c,a,b,j,i,k) = -t3(a,b,c,i,j,k);
                            t3(c,a,b,k,j,i) = -t3(a,b,c,i,j,k);
                        end
                    end
                end
            end
        end
    end
       
    
end
    


