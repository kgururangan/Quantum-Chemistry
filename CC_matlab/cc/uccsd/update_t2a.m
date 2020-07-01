function [t2a] = update_t2a(t1a, t1b, t2a, t2b, t2c, chi1A, chi1B, chi2A, chi2B, chi2C, sys)

    X2a_abij = sys.vA_vvoo... % 1
               +einsum_kg(chi1A.ae,t2a,'ae,ebij->abij')... % 8
               -einsum_kg(chi1A.ae,t2a,'be,eaij->abij')... % 8
               -einsum_kg(chi1A.mi,t2a,'mi,abmj->abij')... % 8
               +einsum_kg(chi1A.mi,t2a,'mj,abmi->abij')... % 8
               +einsum_kg(chi2A.amie_bar,t2a,'amie,ebmj->abij')... % 6
               -einsum_kg(chi2A.amie_bar,t2a,'amje,ebmi->abij')... % 6
               -einsum_kg(chi2A.amie_bar,t2a,'bmie,eamj->abij')... % 6
               +einsum_kg(chi2A.amie_bar,t2a,'bmje,eami->abij')... % 6
               +einsum_kg(chi2B.amie_bar,t2b,'amie,bejm->abij')... % 6
               -einsum_kg(chi2B.amie_bar,t2b,'amje,beim->abij')... % 6
               -einsum_kg(chi2B.amie_bar,t2b,'bmie,aejm->abij')... % 6
               +einsum_kg(chi2B.amie_bar,t2b,'bmje,aeim->abij')... % 6
               +0.5*einsum_kg(chi2A.mnij,t2a,'mnij,abmn->abij')... % 5
               +0.5*einsum_kg(chi2A.abef,t2a,'abef,efij->abij')... % 5
               -einsum_kg(chi2A.mbij,t1a,'mbij,am->abij')... % 8
               +einsum_kg(chi2A.mbij,t1a,'maij,bm->abij')... % 8
               +einsum_kg(chi2A.abej,t1a,'abej,ei->abij')... % 2
               -einsum_kg(chi2A.abej,t1a,'abei,ej->abij'); % 2

% used with chi2A.amie_bar0 and chi2B.amie_bar0... shows the VT2^2 terms
% are not the problem!
% 
%     VAT2A = einsum_kg(sys.vA_oovv,t2a,'mnef,aeim->anif'); % A_{ab}
%     VBT2B = einsum_kg(sys.vB_oovv,t2b,'nmfe,aeim->anif'); % A_{ij}A_{ab}
%     VCT2B = einsum_kg(sys.vC_oovv,t2b,'mnef,aeim->anif'); % A_{ab}

%     X2a_abij = X2a_abij  ...
%                + einsum_kg(VAT2A,t2a,'anif,bfjn->abij') ...
%                - einsum_kg(VAT2A,t2a,'bnif,afjn->abij') ...
%                + einsum_kg(VBT2B,t2a,'anif,bfjn->abij') ...
%                - einsum_kg(VBT2B,t2a,'bnif,afjn->abij') ...
%                - einsum_kg(VBT2B,t2a,'anjf,bfin->abij') ...
%                + einsum_kg(VBT2B,t2a,'bnjf,afin->abij') ...
%                + einsum_kg(VCT2B,t2b,'anif,bfjn->abij') ...
%                - einsum_kg(VCT2B,t2b,'bnif,afjn->abij');

    
    omega = 1;
    for a = 1:sys.Nvir_alpha
        for b = a+1:sys.Nvir_alpha
            for i = 1:sys.Nocc_alpha
                for j = i+1:sys.Nocc_alpha
                    temp = X2a_abij(a,b,i,j)/...
                                       (sys.fa_oo(i,i)+sys.fa_oo(j,j)-sys.fa_vv(a,a)-sys.fa_vv(b,b));
                    t2a(a,b,i,j) = (1-omega)*t2a(a,b,i,j) + omega*temp;                
                    t2a(b,a,i,j) = -t2a(a,b,i,j);
                    t2a(a,b,j,i) = -t2a(a,b,i,j);
                    t2a(b,a,j,i) = t2a(a,b,i,j);
                end
            end
        end
    end


end

