function [t2a] = update_t2a_debug(t1a, t1b, t2a, t2b, t2c, sys, shift)

    D1 = sys.vA_vvoo;
    
    D2 = einsum_kg(sys.vA_vvvo,t1a,'abej,ei->abij');
    D2 = D2 - permute(D2,[1,2,4,3]);
    
    D3 = -einsum_kg(sys.vA_vooo,t1a,'amij,bm->abij');
    D3 = D3 - permute(D3,[2,1,3,4]);
    
    D4 = einsum_kg(einsum_kg(sys.vA_vvvv,t1a,'abef,ei->abif'),t1a,'abif,fj->abij');
    
    D5 = -einsum_kg(einsum_kg(sys.vA_ovvo,t1a,'mbej,ei->mbij'),t1a,'mbij,am->abij');
    D5 = D5 - permute(D5,[2,1,3,4]) - permute(D5,[1,2,4,3]) + permute(D5,[2,1,4,3]);
    
    D6 = einsum_kg(einsum_kg(sys.vA_oooo,t1a,'mnij,bn->mbij'),t1a,'mbij,am->abij');
    
    D7 = einsum_kg(einsum_kg(einsum_kg(sys.vA_ooov,t1a,'mnie,am->anie'),t1a,'anie,bn->abie'),t1a,'abie,ej->abij');
    D7 = D7 - permute(D7,[1,2,4,3]);
    
    D8 = -einsum_kg(einsum_kg(einsum_kg(sys.vA_ovvv,t1a,'mbef,ei->mbif'),t1a,'mbif,fj->mbij'),t1a,'mbij,am->abij');
    D8 = D8 - permute(D8,[2,1,3,4]);
    
    D9 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ei->mnif'),t1a,'mnif,fj->mnij'),t1a,'mnij,am->anij'),t1a,'anij,bn->abij');
    
    MM12 = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9;
    
    % CCS HBar elements
    [HBar_t] = build_hbar_ccs_debug(t1a,t1b,sys);
    H1A = HBar_t.H1A; H1B = HBar_t.H1B; 
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    
    % <ijab|[H(CCS)*T2]_C|0>
    D10 = einsum_kg(H1A.vv-diag(diag(sys.fa_vv)),t2a,'ae,ebij->abij');
    D10 = D10 - permute(D10,[2,1,3,4]);
    
    D11 = -einsum_kg(H1A.oo-diag(diag(sys.fa_oo)),t2a,'mi,abmj->abij');
    D11 = D11 - permute(D11,[1,2,4,3]);
    
    D12 = einsum_kg(H2A.voov,t2a,'amie,ebmj->abij');
    D12 = D12 - permute(D12,[2,1,3,4]) - permute(D12,[1,2,4,3]) + permute(D12,[2,1,4,3]);
    
    D13 = einsum_kg(H2B.voov,t2b,'amie,bejm->abij');
    D13 = D13 - permute(D13,[2,1,3,4]) - permute(D13,[1,2,4,3]) + permute(D13,[2,1,4,3]);
    
    D14 = 0.5*einsum_kg(H2A.vvvv,t2a,'abef,efij->abij');
    
    D15 = 0.5*einsum_kg(H2A.oooo,t2a,'mnij,abmn->abij');
    
    CCS_T2 = D10 + D11 + D12 + D13 + D14 + D15;
    
    % <ijab|[H(CCS)*T2^2]_C|0>
    D16 = einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnef,afin->amie'),t2a,'amie,ebmj->abij');
    D16 = D16 - permute(D16,[1,2,4,3]);
    
    D17 = einsum_kg(einsum_kg(sys.vB_oovv,t2b,'mnef,afin->amie'),t2a,'amie,ebmj->abij');
    D17 = D17 - permute(D17,[2,1,3,4]) - permute(D17,[1,2,4,3]) + permute(D17,[2,1,4,3]);
    
    D18 = 0.25*einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnef,efij->mnij'),t2a,'mnij,abmn->abij');
    
    D19 = -0.5*einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnef,afmn->ae'),t2a,'ae,ebij->abij');
    D19 = D19 - permute(D19,[2,1,3,4]);
    
    D20 = einsum_kg(einsum_kg(sys.vC_oovv,t2b,'mnef,afin->amie'),t2b,'amie,bejm->abij');
    D20 = D20 - permute(D20,[1,2,4,3]);
    
    D21 = -0.5*einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi'),t2a,'mi,abmj->abij');
    D21 = D21 - permute(D21,[1,2,4,3]);
    
    D22 = -einsum_kg(einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ae'),t2a,'ae,ebij->abij');
    D22 = D22 - permute(D22,[2,1,3,4]);
    
    D23 = -einsum_kg(einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi'),t2a,'mi,abmj->abij');
    D23 = D23 - permute(D23,[1,2,4,3]);
    
    CCS_T22 =  D16 + D17 + D18 + D19 + D20 + D21 + D22 + D23;
    
    X2A = MM12 + CCS_T2 + CCS_T22;
    
    
    t2a = X2A;
    
    
%     for i = 1:sys.Nocc_alpha
%         for j = i+1:sys.Nocc_alpha
%             for a = 1:sys.Nvir_alpha
%                 for b = a+1:sys.Nvir_alpha
%                     t2a(a,b,i,j) = X2A(a,b,i,j)/(sys.fa_oo(i,i)+sys.fa_oo(j,j)-sys.fa_vv(a,a)-sys.fa_vv(b,b)-shift);                
%                     t2a(b,a,i,j) = -t2a(a,b,i,j);
%                     t2a(a,b,j,i) = -t2a(a,b,i,j);
%                     t2a(b,a,j,i) = t2a(a,b,i,j);
%                 end
%             end
%         end
%     end
    
end

