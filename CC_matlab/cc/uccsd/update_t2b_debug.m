function [t2b] = update_t2b_debug(t1a, t1b, t2a, t2b, t2c, sys, shift)

    % MM12 contribution
    d1 = sys.vB_vvoo;
    d2 = -einsum_kg(sys.vB_ovoo,t1a,'mbij,am->abij');
    d3 = -einsum_kg(sys.vB_vooo,t1b,'amij,bm->abij');
    d4 = einsum_kg(sys.vB_vvvo,t1a,'abej,ei->abij');
    d5 = einsum_kg(sys.vB_vvov,t1b,'abie,ej->abij');
    d6 = -einsum_kg(einsum_kg(sys.vB_voov,t1b,'amie,ej->amij'),t1b,'amij,bm->abij');
    d7 = -einsum_kg(einsum_kg(sys.vB_ovvo,t1a,'mbej,ei->mbij'),t1a,'mbij,am->abij');
    d8 = einsum_kg(einsum_kg(sys.vB_vvvv,t1b,'abef,fj->abej'),t1a,'abej,ei->abij');
    d9 = einsum_kg(einsum_kg(sys.vB_oooo,t1b,'mnij,bn->mbij'),t1a,'mbij,am->abij');
    d10 = einsum_kg(einsum_kg(einsum_kg(sys.vB_oovo,t1a,'mnej,ei->mnij'),t1b,'mnij,bn->mbij'),t1a,'mbij,am->abij');
    d11 = einsum_kg(einsum_kg(einsum_kg(sys.vB_ooov,t1a,'mnif,am->anif'),t1b,'anif,fj->anij'),t1b,'anij,bn->abij');
    d12 = -einsum_kg(einsum_kg(einsum_kg(sys.vB_ovvv,t1b,'mbef,fj->mbej'),t1a,'mbej,ei->mbij'),t1a,'mbij,am->abij');
    d13 = -einsum_kg(einsum_kg(einsum_kg(sys.vB_vovv,t1b,'anef,fj->anej'),t1a,'anej,ei->anij'),t1b,'anij,bn->abij');
    d14 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ei->mnif'),t1a,'mnif,am->anif'),t1b,'anif,fj->anij'),t1b,'anij,bn->abij');
    d15 = -einsum_kg(einsum_kg(sys.vB_ovov,t1a,'mbif,am->abif'),t1b,'abif,fj->abij');
    d16 = -einsum_kg(einsum_kg(sys.vB_vovo,t1b,'anej,bn->abej'),t1a,'abej,ei->abij');
    
    MM12B = d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8...
            + d9 + d10 + d11 + d12 + d13 + d14 + d15 + d16;

    % CCS HBar elements
    h1A_ov = sys.fa_ov...
             +einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me')...
             +einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');  
         
    h1B_ov = sys.fb_ov...
             +einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me')...
             +einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');

    h1A_oo = sys.fa_oo_masked...
             +einsum_kg(h1A_ov,t1a,'me,ei->mi')...
             +einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi')...
             +einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi');
         
    h1B_oo = sys.fb_oo_masked...
             +einsum_kg(h1B_ov,t1b,'me,ei->mi')...
             +einsum_kg(sys.vB_oovo,t1a,'nmfi,fn->mi')...
             +einsum_kg(sys.vC_ooov,t1b,'mnif,fn->mi');
         
    h1A_vv = sys.fa_vv_masked...
             -einsum_kg(h1A_ov,t1a,'me,am->ae')...
             +einsum_kg(sys.vA_vovv,t1a,'amef,fm->ae')...
             +einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae');
         
    h1B_vv = sys.fb_vv_masked...
             -einsum_kg(h1B_ov,t1b,'me,am->ae')...
             +einsum_kg(sys.vB_ovvv,t1a,'nafe,fn->ae')...
             +einsum_kg(sys.vC_vovv,t1b,'anef,fn->ae');
        
    h2B_oooo =  sys.vB_oooo...
                +einsum_kg(sys.vB_oovo,t1a,'mnej,ei->mnij')...
                +einsum_kg(sys.vB_ooov,t1b,'mnif,fj->mnij')...
                +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ei->mnif'),t1b,'mnif,fj->mnij'); 
            
    h2B_vvvv =  sys.vB_vvvv...
                -einsum_kg(sys.vB_ovvv,t1a,'mbef,am->abef')...
                -einsum_kg(sys.vB_vovv,t1b,'anef,bn->abef')...
                +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,am->anef'),t1b,'anef,bn->abef');
        
    h2A_voov = sys.vA_voov ...
               -einsum_kg(sys.vA_ooov,t1a,'nmie,an->amie')...
               +einsum_kg(sys.vA_vovv,t1a,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie');
            
    h2B_voov = sys.vB_voov ...
               -einsum_kg(sys.vB_ooov,t1a,'nmie,an->amie')...
               +einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie');
           
    h2B_ovov = sys.vB_ovov ...
               +einsum_kg(sys.vB_ovvv,t1a,'mafe,fi->maie')...
               -einsum_kg(sys.vB_ooov,t1b,'mnie,an->maie')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnfe,an->mafe'),t1a,'mafe,fi->maie');
           
    h2B_vovo = sys.vB_vovo ...
               -einsum_kg(sys.vB_oovo,t1a,'nmei,an->amei')...
               +einsum_kg(sys.vB_vovv,t1b,'amef,fi->amei')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmef,fi->nmei'),t1a,'nmei,an->amei');
           
    h2B_ovvo = sys.vB_ovvo ...
               +einsum_kg(sys.vB_ovvv,t1b,'maef,fi->maei')...
               -einsum_kg(sys.vB_oovo,t1b,'mnei,an->maei')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->maei');
           
    h2C_voov = sys.vC_voov ...
               -einsum_kg(sys.vC_oovo,t1b,'mnei,an->amie')...
               +einsum_kg(sys.vC_vovv,t1b,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->amie');
    
    % <ijab|[H(CCS)*T2]_C|0>
    
    D5 = einsum_kg(h1A_vv,t2b,'ae,ebij->abij');
    
    D6 = einsum_kg(h1B_vv,t2b,'be,aeij->abij');
    
    D7 = -einsum_kg(h1A_oo,t2b,'mi,abmj->abij');
    
    D8 = -einsum_kg(h1B_oo,t2b,'mj,abim->abij');
    
    D9 = einsum_kg(h2A_voov,t2b,'amie,ebmj->abij');
    
    D10 = einsum_kg(h2B_voov,t2c,'amie,ebmj->abij');
    
    D11 = einsum_kg(h2B_ovvo,t2a,'mbej,aeim->abij');
    
    D12 = einsum_kg(h2C_voov,t2b,'bmje,aeim->abij');
    
    D13 = -einsum_kg(h2B_ovov,t2b,'mbie,aemj->abij');
    
    D14 = -einsum_kg(h2B_vovo,t2b,'amej,ebim->abij');
    
    D15 = einsum_kg(h2B_oooo,t2b,'mnij,abmn->abij');
    
    D16 = einsum_kg(h2B_vvvv,t2b,'abef,efij->abij');
    
    CCS_T2 = D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14 + D15 + D16;
    
    % <ijab|[H(CCS)*T2^2]_C|0>
    
    D17 = einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnef,aeim->anif'),t2b,'anif,fbnj->abij');
    
    D18 = einsum_kg(einsum_kg(sys.vB_oovv,t2b,'nmfe,aeim->anif'),t2b,'anif,fbnj->abij');
    
    D19 = einsum_kg(einsum_kg(sys.vB_oovv,t2a,'mnef,aeim->anif'),t2c,'anif,fbnj->abij');
    
    D20 = einsum_kg(einsum_kg(sys.vC_oovv,t2b,'mnef,aeim->anif'),t2c,'anif,fbnj->abij');
    
    D21 = einsum_kg(einsum_kg(sys.vB_oovv,t2b,'mnef,afmj->anej'),t2b,'anej,ebin->abij');
    
    D22 = -0.5*einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnef,afmn->ae'),t2b,'ae,ebij->abij');
    
    D23 = -einsum_kg(einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ae'),t2b,'ae,ebij->abij');
    
    D24 = -einsum_kg(einsum_kg(sys.vB_oovv,t2b,'nmfe,fbnm->be'),t2b,'be,aeij->abij');
    
    D25 = -0.5*einsum_kg(einsum_kg(sys.vC_oovv,t2c,'mnef,fbnm->be'),t2b,'be,aeij->abij');
    
    D26 = -0.5*einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi'),t2b,'mi,abmj->abij');
    
    D27 = -einsum_kg(einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi'),t2b,'mi,abmj->abij');
    
    D28 = -einsum_kg(einsum_kg(sys.vB_oovv,t2b,'nmfe,fenj->mj'),t2b,'mj,abim->abij');
    
    D29 = -0.5*einsum_kg(einsum_kg(sys.vC_oovv,t2c,'mnef,efjn->mj'),t2b,'mj,abim->abij');
    
    D30 = einsum_kg(einsum_kg(sys.vB_oovv,t2b,'mnef,efij->mnij'),t2b,'mnij,abmn->abij');
     
    CCS_T22 = D17 + D18 + D19 + D20 + D21 + D22 + D23 + D24 + D25 + D26 + D27 + D28 + D29 + D30;
    
    X2B = MM12B + CCS_T2 + CCS_T22;
    
    for i = 1:sys.Nocc_alpha
        for j = 1:sys.Nocc_beta
            for a = 1:sys.Nvir_alpha
                for b = 1:sys.Nvir_beta
                    denom = (sys.fa_oo(i,i)+sys.fb_oo(j,j)-sys.fa_vv(a,a)-sys.fb_vv(b,b)-shift); 
                    t2b(a,b,i,j) = X2B(a,b,i,j) / denom;
                end
            end
        end
    end
    
end
