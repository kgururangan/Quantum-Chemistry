function [MM23] = build_MM23_debug(t1,t2,sys)

    % very, very wrong at the moment...

    fprintf('MMCC(2,3) debug construction... ')
    
    tic

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
  
                                               
%     d1 = sys.Vvooo;
%     f1 = sys.Vvooo;
%     get_error(d1,f1)
%     
%     d2 = +einsum_kg(sys.Vvovo,t1,'amej,ei->amij') - einsum_kg(sys.Vvovo,t1,'amei,ej->amij');
%     f2 = +einsum_kg(sys.Vvoov,t1,'amie,ej->amij') - einsum_kg(sys.Vvoov,t1,'amje,ei->amij');
%     get_error(d2,f2)
%     
%     d3 = -einsum_kg(sys.Voooo,t1,'nmij,an->amij');
%     f3 = -einsum_kg(sys.Voooo,t1,'nmij,an->amij');
%     get_error(d3,f3)
%     
%     d4 = +einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fi->amie'),t1,'amie,ej->amij');
%     f4 = +einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fi->amie'),t1,'amie,ej->amij');
%     get_error(d4,f4)
%     
%     d5 = -einsum_kg(einsum_kg(sys.Voovo,t1,'nmej,an->amej'),t1,'amej,ei->amij') + einsum_kg(einsum_kg(sys.Voovo,t1,'nmei,an->amei'),t1,'amei,ej->amij');
%     f5 = -einsum_kg(einsum_kg(sys.Vooov,t1,'nmie,ej->nmij'),t1,'nmij,an->amij') + einsum_kg(einsum_kg(sys.Vooov,t1,'nmje,ei->nmij'),t1,'nmij,an->amij');
%     get_error(d5,f5)
%     
%     d6 = 
%     f6 = -einsum_kg(einsum_kg(sys.Voovo,t1,'nmfj,an->amfj'),t1,'amfj,fi->amij') + einsum_kg(einsum_kg(sys.Voovo,t1,'nmfi,an->amfi'),t1,'amfi,fj->amij');
% 
%     
% %     H2{2,1,1,1} = sys.Vvooo ...
% %                  + einsum_kg(sys.Vvovo,t1,'amej,ei->amij') - einsum_kg(sys.Vvovo,t1,'amei,ej->amij') ...
% %                  - einsum_kg(sys.Voooo,t1,'nmij,an->amij') ...
% %                  + einsum_kg(sys.fov,t2,'me,aeij->amij') ...
% %                  + 0.5*einsum_kg(sys.Vvovv,t2,'amfe,feij->amij') ...
% %                  + einsum_kg(sys.Voovo,t2,'nmej,aein->amij') - einsum_kg(sys.Voovo,t2,'nmei,aejn->amij') ...
% %                  - einsum_kg(einsum_kg(sys.Voovo,t1,'nmej,an->amej'),t1,'amej,ei->amij') + einsum_kg(einsum_kg(sys.Voovo,t1,'nmei,an->amei'),t1,'amei,ej->amij') ...
% %                  + einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fi->amie'),t1,'amie,ej->amij') ...
% %                  - 0.5*einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,an->amef'),t2,'amef,feij->amij') ...
% %                  + einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,ej->nmfj'),t2,'nmfj,afin->amij') - einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,ei->nmfi'),t2,'nmfi,afjn->amij') ...
% %                  + einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,fn->me'),t2,'me,aeij->amij') ...
% %                  - einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,an->amfe'),t1,'amfe,fi->amie'),t1,'amie,ej->amij'); 
% %                                                
%                                                
%     h_vooo = -(-sys.Vvooo...
%                -einsum_kg(sys.Vvoov,t1,'amie,ej->amij')+einsum_kg(sys.Vvoov,t1,'amje,ei->amij')...
%                +einsum_kg(sys.Voooo,t1,'nmij,an->amij')...
%                -einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fi->amie'),t1,'amie,ej->amij')...
%                +einsum_kg(einsum_kg(sys.Vooov,t1,'nmie,ej->nmij'),t1,'nmij,an->amij') - einsum_kg(einsum_kg(sys.Vooov,t1,'nmje,ei->nmij'),t1,'nmij,an->amij')...
%                +einsum_kg(einsum_kg(sys.Voovo,t1,'nmfj,an->amfj'),t1,'amfj,fi->amij') - einsum_kg(einsum_kg(sys.Voovo,t1,'nmfi,an->amfi'),t1,'amfi,fj->amij')...
%                +einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,ej->nmfj'),t1,'nmfj,an->amfj'),t1,'amfj,fi->amij')...
%                -0.5*einsum_kg(sys.Vvovv,t2,'amfe,feij->amij')...
%                -einsum_kg(sys.Vooov,t2,'mnjf,afin->amij') + einsum_kg(sys.Vooov,t2,'mnif,afjn->amij')...
%                +einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'nmfe,fi->nmie'),t1,'nmie,an->amie'),t1,'amie,ej->amij')...
%                -einsum_kg(sys.fov,t2,'me,aeij->amij')...
%                -einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fn->me'),t2,'me,aeij->amij'));
%     
%     h_vvvo = sys.Vvvvo...
%              -einsum_kg(sys.Vovvo,t1,'mcek,bm->bcek') + einsum_kg(sys.Vovvo,t1,'mbek,cm->bcek')...
%              +einsum_kg(sys.Vvvvv,t1,'bcef,fk->bcek')...
%              +einsum_kg(einsum_kg(sys.Voovo,t1,'mnek,cn->mcek'),t1,'mcek,bm->bcek')...
%              -einsum_kg(einsum_kg(sys.Vovvv,t1,'mcef,fk->mcek'),t1,'mcek,bm->bcek') + einsum_kg(einsum_kg(sys.Vovvv,t1,'mbef,fk->mbek'),t1,'mbek,cm->bcek')...
%              -einsum_kg(einsum_kg(sys.Vvovv,t1,'bnef,cn->bcef'),t1,'bcef,fk->bcek') + einsum_kg(einsum_kg(sys.Vvovv,t1,'cnef,bn->bcef'),t1,'bcef,fk->bcek')...
%              +einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,cn->mcef'),t1,'mcef,bm->bcef'),t1,'bcef,fk->bcek')...
%              +0.5*einsum_kg(sys.Voovo,t2,'mnek,bcmn->bcek')...
%              +einsum_kg(sys.Vvovv,t2,'bnef,fcnk->bcek') - einsum_kg(sys.Vvovv,t2,'cnef,fbnk->bcek')...
%              +einsum_kg(einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,bm->bnef'),t1,'bnef,cn->bcef'),t1,'bcef,fk->bcek')...
%              -einsum_kg(sys.fov,t2,'me,bcmk->bcek')...
%              -einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fn->me'),t2,'me,bcmk->bcek');
%     h_vvov = -permute(h_vvvo,[1,2,4,3]);
%     

%        

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
    
%     MM23 =   D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 ...
%           + D13 + D14 + D15 + D16 + D17 + D18 + D19 + D22 + D23 + D24 + D25 + D26;
% 
%     % A(k/ij)A(a/bc)
%     MM23 = MM23 - permute(MM23,[2,1,3,4,5,6]) - permute(MM23,[3,2,1,4,5,6]);
%     MM23 = MM23 - permute(MM23,[1,2,3,6,5,4]) - permute(MM23,[1,2,3,4,6,5]);

    fprintf('finished in %4.2f s\n',toc)
end
