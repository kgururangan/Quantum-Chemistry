function [t2b] = update_t2b(t1a, t1b, t2a, t2b, t2c, chi1A, chi1B, chi2A, chi2B, chi2C, sys, shift)

    X2b_abij = sys.vB_vvoo ... % 1
               -einsum_kg(chi1A.mi,t2b,'mi,abmj->abij') ... % 8
               -einsum_kg(chi1B.mj,t2b,'mj,abim->abij')... % 8
               +einsum_kg(chi1A.ae,t2b,'ae,ebij->abij')... % 8
               +einsum_kg(chi1B.be,t2b,'be,aeij->abij')... % 8
               +einsum_kg(chi2A.amie,t2b,'amie,ebmj->abij')... % 6
               +einsum_kg(chi2B.amie,t2c,'amie,ebmj->abij')... % 6
               +einsum_kg(chi2B.mbej,t2a,'mbej,aeim->abij')... % 4
               +einsum_kg(chi2C.bmje,t2b,'bmje,aeim->abij')... % 4
               +einsum_kg(chi2B.mnij,t2b,'mnij,abmn->abij')... % 5
               +einsum_kg(chi2B.abef,t2b,'abef,efij->abij')... % 4
               -einsum_kg(chi2B.amej,t2b,'amej,ebim->abij')... % 4
               -einsum_kg(chi2B.mbie,t2b,'mbie,aemj->abij')... % 5
               +einsum_kg(chi2B.abej,t1a,'abej,ei->abij')... % 7
               +einsum_kg(chi2B.abie,t1b,'abie,ej->abij')... % 5
               -einsum_kg(chi2B.mbij,t1a,'mbij,am->abij')... % 2
               -einsum_kg(sys.vB_vooo,t1b,'amij,bm->abij'); % 1


%      VAT2A = einsum_kg(sys.vA_oovv,t2a,'mnef,afin->amie'); 
%      VBT2B = einsum_kg(sys.vB_oovv,t2b,'mnef,afin->amie'); 
%      VCT2B = einsum_kg(sys.vC_oovv,t2b,'nmfe,afin->amie'); 
%      VBT2A = einsum_kg(sys.vB_oovv,t2a,'nmfe,afin->amie');
% % 
%     X2b_abij = X2b_abij  ...
%               +einsum_kg(VAT2A,t2b,'amie,ebmj->abij') ...
%               +einsum_kg(VBT2B,t2b,'amie,ebmj->abij') ...
%               +einsum_kg(VCT2B+VBT2A,t2c,'amie,ebmj->abij');
              
           
          
           
     omega = 1;
     for a = 1:sys.Nvir_alpha
         for b = 1:sys.Nvir_beta
             for i = 1:sys.Nocc_alpha
                 for j = 1:sys.Nocc_beta
                     
                     temp =     X2b_abij(a,b,i,j)/...
                                        (sys.fa_oo(i,i)+sys.fb_oo(j,j)-sys.fa_vv(a,a)-sys.fb_vv(b,b)-shift);  
                     t2b(a,b,i,j) = (1-omega)*t2b(a,b,i,j) + omega*temp;
                 end
             end
         end
     end




end

