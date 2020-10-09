function [t2b] = update_t2b_ccsdt3(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift)

    % MM12 contribution
    d1 = sys.vB_ovoo;
    d2 = einsum_kg(sys.vB_ovvo,t1a,'mbej,ei->mbij');
    d3 = -einsum_kg(sys.vB_oooo,t1b,'mnij,bn->mbij');
    d4 = -einsum_kg(einsum_kg(sys.vB_ooov,t1b,'mnif,bn->mbif'),t1b,'mbif,fj->mbij');
    d5 = -einsum_kg(einsum_kg(sys.vB_oovo,t1b,'mnej,bn->mbej'),t1a,'mbej,ei->mbij');
    d6 = einsum_kg(einsum_kg(sys.vB_ovvv,t1b,'mbef,fj->mbej'),t1a,'mbej,ei->mbij');
    
    h2B_ovoo = d1 + d2 + d3 + d4 + d5 + d6;
    
    d1 = sys.vB_vooo;
    d2 = einsum_kg(sys.vB_voov,t1b,'amif,fj->amij');
    d3 = -einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmef,an->amef'),t1a,'amef,ei->amif'),t1b,'amif,fj->amij');
    d4 = einsum_kg(einsum_kg(sys.vB_vovv,t1b,'amef,fj->amej'),t1a,'amej,ei->amij');
    
    h2B_vooo = d1 + d2 + d3 + d4;
    
    d1 = sys.vB_vvvo;
    d2 = einsum_kg(sys.vB_vvvv,t1b,'abef,fj->abej');
    d3 = -einsum_kg(sys.vB_vovo,t1b,'anej,bn->abej');
    
    h2B_vvvo = d1 + d2 + d3;
     
    d1 = sys.vB_vvov;
    d2 = -einsum_kg(sys.vB_ovov,t1a,'mbie,am->abie');
    
    h2B_vvov = d1 + d2;
     
    D1 = -einsum_kg(h2B_ovoo,t1a,'mbij,am->abij');
    D2 = -einsum_kg(h2B_vooo,t1b,'amij,bm->abij');
    D3 = einsum_kg(h2B_vvvo,t1a,'abej,ei->abij');
    D4 = einsum_kg(h2B_vvov,t1b,'abie,ej->abij');
    
    MM12B = sys.vB_vvoo + D1 + D2 + D3 + D4;

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
    
    % <ijab|[H(CCS)*T2]_C + [H(CCS)*T2^2]_C|0>
    VT2_1A_vv = -0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,afmn->ae') - einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ae');
    VT2_1B_vv = -einsum_kg(sys.vB_oovv,t2b,'nmfe,fbnm->be') - 0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,fbnm->be');
    VT2_1A_oo = 0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi') + einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi');
    VT2_1B_oo = einsum_kg(sys.vB_oovv,t2b,'nmfe,fenj->mj') + 0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,efjn->mj');
    
    VT2_2A_voov = einsum_kg(sys.vA_oovv,t2a,'mnef,aeim->anif') + einsum_kg(sys.vB_oovv,t2b,'nmfe,aeim->anif');
    VT2_2B_voov = einsum_kg(sys.vB_oovv,t2a,'mnef,aeim->anif') + einsum_kg(sys.vC_oovv,t2b,'mnef,aeim->anif');
    VT2_2B_oooo = einsum_kg(sys.vB_oovv,t2b,'mnef,efij->mnij');
    VT2_2B_vovo = -einsum_kg(sys.vB_oovv,t2b,'mnef,afmj->anej');
    
    D5 = einsum_kg(h1A_vv + VT2_1A_vv,t2b,'ae,ebij->abij');
    
    D6 = einsum_kg(h1B_vv + VT2_1B_vv,t2b,'be,aeij->abij');
    
    D7 = -einsum_kg(h1A_oo + VT2_1A_oo,t2b,'mi,abmj->abij');
    
    D8 = -einsum_kg(h1B_oo + VT2_1B_oo,t2b,'mj,abim->abij');
    
    D9 = einsum_kg(h2A_voov + VT2_2A_voov,t2b,'amie,ebmj->abij');
    
    D10 = einsum_kg(h2B_voov + VT2_2B_voov,t2c,'amie,ebmj->abij');
    
    D11 = einsum_kg(h2B_ovvo,t2a,'mbej,aeim->abij');
    
    D12 = einsum_kg(h2C_voov,t2b,'bmje,aeim->abij');
    
    D13 = -einsum_kg(h2B_ovov,t2b,'mbie,aemj->abij');
    
    D14 = -einsum_kg(h2B_vovo + VT2_2B_vovo,t2b,'amej,ebim->abij');
    
    D15 = einsum_kg(h2B_oooo + VT2_2B_oooo,t2b,'mnij,abmn->abij');
    
    D16 = einsum_kg(h2B_vvvv,t2b,'abef,efij->abij');
    
    CCS_T2 = D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14 + D15 + D16;
    
    X2B = MM12B + CCS_T2;
    
    t2b = zeros(sys.Nvir_alpha,sys.Nvir_beta,sys.Nocc_alpha,sys.Nocc_beta);
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
    
    
    % CCSDT part
    PA = sys.PA; PB = sys.PB; HA = sys.HA; HB = sys.HB;
    pA = sys.pA; pB = sys.pB; hA = sys.hA; hB = sys.hB;
    
    h2A_ooov = sys.vA_ooov + einsum_kg(sys.vA_oovv,t1a,'mnfe,fi->mnie');
    h2A_vovv = sys.vA_vovv - einsum_kg(sys.vA_oovv,t1a,'mnfe,an->amef');
    h2B_vovv = sys.vB_vovv - einsum_kg(sys.vB_oovv,t1a,'nmef,an->amef');
    h2B_ovvv = sys.vB_ovvv - einsum_kg(sys.vB_oovv,t1b,'mnfe,an->mafe');
    h2B_ooov = sys.vB_ooov + einsum_kg(sys.vB_oovv,t1a,'mnfe,fi->mnie');
    h2B_oovo = sys.vB_oovo + einsum_kg(sys.vB_oovv,t1b,'nmef,fi->nmei');
    h2C_ooov = sys.vC_ooov + einsum_kg(sys.vC_oovv,t1b,'mnfe,fi->mnie');
    h2C_vovv = sys.vC_vovv - einsum_kg(sys.vC_oovv,t1b,'mnfe,an->amef');
    
    % IJAB
    D1 = -0.5*einsum_kg(h2A_ooov(HA,HA,HA,PA),t3b,'mnif,afbmnj->abij');
    D2 = -einsum_kg(h2B_oovo(HA,HB,PA,HB),t3b,'nmfj,afbinm->abij');
    D3 = -0.5*einsum_kg(h2C_ooov(HB,HB,HB,PB),t3c,'mnjf,afbinm->abij');
    D4 = -einsum_kg(h2B_ooov(HA,HB,HA,PB),t3c,'mnif,afbmnj->abij');
    D5 = +0.5*einsum_kg(h2A_vovv(PA,HA,PA,PA),t3b,'anef,efbinj->abij');
    D6 = +einsum_kg(h2B_vovv(PA,HB,PA,PB),t3c,'anef,efbinj->abij');
    D7 = +einsum_kg(h2B_ovvv(HA,PB,PA,PB),t3b,'nbfe,afeinj->abij');
    D8 = +0.5*einsum_kg(h2C_vovv(PB,HB,PB,PB),t3c,'bnef,afeinj->abij');
    D9 = +einsum_kg(h1A_ov(HA,PA),t3b,'me,aebimj->abij');
    D10 = +einsum_kg(h1B_ov(HB,PB),t3c,'me,aebimj->abij');
    
    X2B_IJAB = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10;
    
    t2b_ABIJ = zeros(sys.Nact_p_alpha,sys.Nact_p_beta,sys.Nact_h_alpha,sys.Nact_h_beta);
    for i = 1:sys.Nact_h_alpha
        for j = 1:sys.Nact_h_beta
            for a = 1:sys.Nact_p_alpha
                for b = 1:sys.Nact_p_beta
                    denom = (sys.fa_HH(i,i)+sys.fb_HH(j,j)-sys.fa_PP(a,a)-sys.fb_PP(b,b)-shift); 
                    t2b_ABIJ(a,b,i,j) = X2B_IJAB(a,b,i,j) / denom;
                end
            end
        end
    end
    t2b(PA,PB,HA,HB) = t2b(PA,PB,HA,HB) + t2b_ABIJ;
    
    % iJAB
    D1 = -0.5*einsum_kg(h2A_ooov(HA,HA,hA,PA),t3b,'mnif,afbmnj->abij');
    D4 = -einsum_kg(h2B_ooov(HA,HB,hA,PB),t3c,'mnif,afbmnj->abij');
    
    X2B_iJAB = D1 + D4;
    
    t2b_ABiJ = zeros(sys.Nact_p_alpha,sys.Nact_p_beta,sys.Nunact_h_alpha,sys.Nact_h_beta);
    for i = 1:sys.Nunact_h_alpha
        for j = 1:sys.Nact_h_beta
            for a = 1:sys.Nact_p_alpha
                for b = 1:sys.Nact_p_beta
                    denom = (sys.fa_hh(i,i)+sys.fb_HH(j,j)-sys.fa_PP(a,a)-sys.fb_PP(b,b)-shift); 
                    t2b_ABiJ(a,b,i,j) = X2B_iJAB(a,b,i,j) / denom;
                end
            end
        end
    end
    t2b(PA,PB,hA,HB) = t2b(PA,PB,hA,HB) + t2b_ABiJ;
    
    % IjAB
    D2 = -einsum_kg(h2B_oovo(HA,HB,PA,hB),t3b,'nmfj,afbinm->abij');
    D3 = -0.5*einsum_kg(h2C_ooov(HB,HB,hB,PB),t3c,'mnjf,afbinm->abij');
    
    X2B_IjAB = D2 + D3;
    
    t2b_ABIj = zeros(sys.Nact_p_alpha,sys.Nact_p_beta,sys.Nact_h_alpha,sys.Nunact_h_beta);
    for i = 1:sys.Nact_h_alpha
        for j = 1:sys.Nunact_h_beta
            for a = 1:sys.Nact_p_alpha
                for b = 1:sys.Nact_p_beta
                    denom = (sys.fa_HH(i,i)+sys.fb_hh(j,j)-sys.fa_PP(a,a)-sys.fb_PP(b,b)-shift); 
                    t2b_ABIj(a,b,i,j) = X2B_IjAB(a,b,i,j) / denom;
                end
            end
        end
    end
    t2b(PA,PB,HA,hB) = t2b(PA,PB,HA,hB) + t2b_ABIj;
    
    % IJaB
    D5 = +0.5*einsum_kg(h2A_vovv(pA,HA,PA,PA),t3b,'anef,efbinj->abij');
    D6 = +einsum_kg(h2B_vovv(pA,HB,PA,PB),t3c,'anef,efbinj->abij');
    
    X2B_IJaB = D5 + D6;
    
    t2b_aBIJ = zeros(sys.Nunact_p_alpha,sys.Nact_p_beta,sys.Nact_h_alpha,sys.Nact_h_beta);
    for i = 1:sys.Nact_h_alpha
        for j = 1:sys.Nact_h_beta
            for a = 1:sys.Nunact_p_alpha
                for b = 1:sys.Nact_p_beta
                    denom = (sys.fa_HH(i,i)+sys.fb_HH(j,j)-sys.fa_pp(a,a)-sys.fb_PP(b,b)-shift); 
                    t2b_aBIJ(a,b,i,j) = X2B_IJaB(a,b,i,j) / denom;
                end
            end
        end
    end
    t2b(pA,PB,HA,HB) = t2b(pA,PB,HA,HB) + t2b_aBIJ;
    
    % IJAb
    D7 = +einsum_kg(h2B_ovvv(HA,pB,PA,PB),t3b,'nbfe,afeinj->abij');
    D8 = +0.5*einsum_kg(h2C_vovv(pB,HB,PB,PB),t3c,'bnef,afeinj->abij');
    
    X2B_IJAb = D7 + D8;
    
    t2b_AbIJ = zeros(sys.Nact_p_alpha,sys.Nunact_p_beta,sys.Nact_h_alpha,sys.Nact_h_beta);
    for i = 1:sys.Nact_h_alpha
        for j = 1:sys.Nact_h_beta
            for a = 1:sys.Nact_p_alpha
                for b = 1:sys.Nunact_p_beta
                    denom = (sys.fa_HH(i,i)+sys.fb_HH(j,j)-sys.fa_PP(a,a)-sys.fb_pp(b,b)-shift); 
                    t2b_AbIJ(a,b,i,j) = X2B_IJAb(a,b,i,j) / denom;
                end
            end
        end
    end
    t2b(PA,pB,HA,HB) = t2b(PA,pB,HA,HB) + t2b_AbIJ;
    
end

