function [X2B] = update_t2b_proj1(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys,shift)

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

    h1A_oo = sys.fa_oo...
             +einsum_kg(h1A_ov,t1a,'me,ei->mi')...
             +einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi')...
             +einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi');
         
    h1B_oo = sys.fb_oo...
             +einsum_kg(h1B_ov,t1b,'me,ei->mi')...
             +einsum_kg(sys.vB_oovo,t1a,'nmfi,fn->mi')...
             +einsum_kg(sys.vC_ooov,t1b,'mnif,fn->mi');
         
    h1A_vv = sys.fa_vv...
             -einsum_kg(h1A_ov,t1a,'me,am->ae')...
             +einsum_kg(sys.vA_vovv,t1a,'amef,fm->ae')...
             +einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae');
         
    h1B_vv = sys.fb_vv...
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
    
    % CCSDT part
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

    H1A.ov = h1A_ov;
    H1B.ov = h1B_ov;
    H2A.ooov = sys.vA_ooov + einsum_kg(sys.vA_oovv,t1a,'mnfe,fi->mnie');
    H2A.vovv = sys.vA_vovv - einsum_kg(sys.vA_oovv,t1a,'mnfe,an->amef');
    H2B.vovv = sys.vB_vovv - einsum_kg(sys.vB_oovv,t1a,'nmef,an->amef');
    H2B.ovvv = sys.vB_ovvv - einsum_kg(sys.vB_oovv,t1b,'mnfe,an->mafe');
    H2B.ooov = sys.vB_ooov + einsum_kg(sys.vB_oovv,t1a,'mnfe,fi->mnie');
    H2B.oovo = sys.vB_oovo + einsum_kg(sys.vB_oovv,t1b,'nmef,fi->nmei');
    H2C.ooov = sys.vC_ooov + einsum_kg(sys.vC_oovv,t1b,'mnfe,fi->mnie');
    H2C.vovv = sys.vC_vovv - einsum_kg(sys.vC_oovv,t1b,'mnfe,an->amef');
    

D1 = -0.5*einsum_kg(H2A.ooov(hA,hA,HA,pA),T3B.PpPhhH,'mnIf,AfBmnJ->ABIJ')...
-0.5*einsum_kg(H2A.ooov(hA,hA,HA,PA),T3B.PPPhhH,'mnIF,AFBmnJ->ABIJ')...
+einsum_kg(H2A.ooov(HA,hA,HA,pA),T3B.PpPhHH,'MnIf,AfBnMJ->ABIJ')...
+einsum_kg(H2A.ooov(HA,hA,HA,PA),T3B.PPPhHH,'MnIF,AFBnMJ->ABIJ')...
-0.5*einsum_kg(H2A.ooov(HA,HA,HA,pA),T3B.PpPHHH,'MNIf,AfBMNJ->ABIJ')...
-0.5*einsum_kg(H2A.ooov(HA,HA,HA,PA),T3B.PPPHHH,'MNIF,AFBMNJ->ABIJ');

D2 = +einsum_kg(H2B.oovo(hA,hB,pA,HB),T3B.PpPhHh,'nmfJ,AfBnIm->ABIJ')...
+einsum_kg(H2B.oovo(hA,hB,PA,HB),T3B.PPPhHh,'nmFJ,AFBnIm->ABIJ')...
+einsum_kg(H2B.oovo(hA,HB,pA,HB),T3B.PpPhHH,'nMfJ,AfBnIM->ABIJ')...
+einsum_kg(H2B.oovo(hA,HB,PA,HB),T3B.PPPhHH,'nMFJ,AFBnIM->ABIJ')...
-einsum_kg(H2B.oovo(HA,hB,pA,HB),T3B.PpPHHh,'NmfJ,AfBINm->ABIJ')...
-einsum_kg(H2B.oovo(HA,hB,PA,HB),T3B.PPPHHh,'NmFJ,AFBINm->ABIJ')...
-einsum_kg(H2B.oovo(HA,HB,pA,HB),T3B.PpPHHH,'NMfJ,AfBINM->ABIJ')...
-einsum_kg(H2B.oovo(HA,HB,PA,HB),T3B.PPPHHH,'NMFJ,AFBINM->ABIJ');

D3 = +0.5*einsum_kg(H2C.ooov(hB,hB,HB,pB),T3C.PPpHhh,'mnJf,ABfInm->ABIJ')...
-0.5*einsum_kg(H2C.ooov(hB,hB,HB,PB),T3C.PPPHhh,'mnJF,AFBInm->ABIJ')...
+einsum_kg(H2C.ooov(HB,hB,HB,pB),T3C.PPpHhH,'MnJf,ABfInM->ABIJ')...
-einsum_kg(H2C.ooov(HB,hB,HB,PB),T3C.PPPHhH,'MnJF,AFBInM->ABIJ')...
+0.5*einsum_kg(H2C.ooov(HB,HB,HB,pB),T3C.PPpHHH,'MNJf,ABfINM->ABIJ')...
-0.5*einsum_kg(H2C.ooov(HB,HB,HB,PB),T3C.PPPHHH,'MNJF,AFBINM->ABIJ');

D4 = +einsum_kg(H2B.ooov(hA,hB,HA,pB),T3C.PPphhH,'mnIf,ABfmnJ->ABIJ')...
-einsum_kg(H2B.ooov(hA,hB,HA,PB),T3C.PPPhhH,'mnIF,AFBmnJ->ABIJ')...
+einsum_kg(H2B.ooov(hA,HB,HA,pB),T3C.PPphHH,'mNIf,ABfmNJ->ABIJ')...
-einsum_kg(H2B.ooov(hA,HB,HA,PB),T3C.PPPhHH,'mNIF,AFBmNJ->ABIJ')...
+einsum_kg(H2B.ooov(HA,hB,HA,pB),T3C.PPpHhH,'MnIf,ABfMnJ->ABIJ')...
-einsum_kg(H2B.ooov(HA,hB,HA,PB),T3C.PPPHhH,'MnIF,AFBMnJ->ABIJ')...
+einsum_kg(H2B.ooov(HA,HB,HA,pB),T3C.PPpHHH,'MNIf,ABfMNJ->ABIJ')...
-einsum_kg(H2B.ooov(HA,HB,HA,PB),T3C.PPPHHH,'MNIF,AFBMNJ->ABIJ');

D5 = -0.5*einsum_kg(H2A.vovv(PA,hA,pA,pA),T3B.ppPhHH,'Anef,efBnIJ->ABIJ')...
-einsum_kg(H2A.vovv(PA,hA,PA,pA),T3B.PpPhHH,'AnEf,EfBnIJ->ABIJ')...
-0.5*einsum_kg(H2A.vovv(PA,hA,PA,PA),T3B.PPPhHH,'AnEF,EFBnIJ->ABIJ')...
+0.5*einsum_kg(H2A.vovv(PA,HA,pA,pA),T3B.ppPHHH,'ANef,efBINJ->ABIJ')...
+einsum_kg(H2A.vovv(PA,HA,PA,pA),T3B.PpPHHH,'ANEf,EfBINJ->ABIJ')...
+0.5*einsum_kg(H2A.vovv(PA,HA,PA,PA),T3B.PPPHHH,'ANEF,EFBINJ->ABIJ');

D6 = -einsum_kg(H2B.vovv(PA,hB,pA,pB),T3C.pPpHhH,'Anef,eBfInJ->ABIJ')...
+einsum_kg(H2B.vovv(PA,hB,pA,PB),T3C.pPPHhH,'AneF,eFBInJ->ABIJ')...
-einsum_kg(H2B.vovv(PA,hB,PA,pB),T3C.PPpHhH,'AnEf,EBfInJ->ABIJ')...
+einsum_kg(H2B.vovv(PA,hB,PA,PB),T3C.PPPHhH,'AnEF,EFBInJ->ABIJ')...
-einsum_kg(H2B.vovv(PA,HB,pA,pB),T3C.pPpHHH,'ANef,eBfINJ->ABIJ')...
+einsum_kg(H2B.vovv(PA,HB,pA,PB),T3C.pPPHHH,'ANeF,eFBINJ->ABIJ')...
-einsum_kg(H2B.vovv(PA,HB,PA,pB),T3C.PPpHHH,'ANEf,EBfINJ->ABIJ')...
+einsum_kg(H2B.vovv(PA,HB,PA,PB),T3C.PPPHHH,'ANEF,EFBINJ->ABIJ');

D7 = -einsum_kg(H2B.ovvv(hA,PB,pA,pB),T3B.PpphHH,'nBfe,AfenIJ->ABIJ')...
-einsum_kg(H2B.ovvv(hA,PB,pA,PB),T3B.PpPhHH,'nBfE,AfEnIJ->ABIJ')...
-einsum_kg(H2B.ovvv(hA,PB,PA,pB),T3B.PPphHH,'nBFe,AFenIJ->ABIJ')...
-einsum_kg(H2B.ovvv(hA,PB,PA,PB),T3B.PPPhHH,'nBFE,AFEnIJ->ABIJ')...
+einsum_kg(H2B.ovvv(HA,PB,pA,pB),T3B.PppHHH,'NBfe,AfeINJ->ABIJ')...
+einsum_kg(H2B.ovvv(HA,PB,pA,PB),T3B.PpPHHH,'NBfE,AfEINJ->ABIJ')...
+einsum_kg(H2B.ovvv(HA,PB,PA,pB),T3B.PPpHHH,'NBFe,AFeINJ->ABIJ')...
+einsum_kg(H2B.ovvv(HA,PB,PA,PB),T3B.PPPHHH,'NBFE,AFEINJ->ABIJ');

D8 = +0.5*einsum_kg(H2C.vovv(PB,hB,pB,pB),T3C.PppHhH,'Bnef,AfeInJ->ABIJ')...
-einsum_kg(H2C.vovv(PB,hB,PB,pB),T3C.PPpHhH,'BnEf,AEfInJ->ABIJ')...
+0.5*einsum_kg(H2C.vovv(PB,hB,PB,PB),T3C.PPPHhH,'BnEF,AFEInJ->ABIJ')...
+0.5*einsum_kg(H2C.vovv(PB,HB,pB,pB),T3C.PppHHH,'BNef,AfeINJ->ABIJ')...
-einsum_kg(H2C.vovv(PB,HB,PB,pB),T3C.PPpHHH,'BNEf,AEfINJ->ABIJ')...
+0.5*einsum_kg(H2C.vovv(PB,HB,PB,PB),T3C.PPPHHH,'BNEF,AFEINJ->ABIJ');

D9 = -einsum_kg(H1A.ov(hA,pA),T3B.PpPhHH,'me,AeBmIJ->ABIJ')...
-einsum_kg(H1A.ov(hA,PA),T3B.PPPhHH,'mE,AEBmIJ->ABIJ')...
+einsum_kg(H1A.ov(HA,pA),T3B.PpPHHH,'Me,AeBIMJ->ABIJ')...
+einsum_kg(H1A.ov(HA,PA),T3B.PPPHHH,'ME,AEBIMJ->ABIJ');

D10 = -einsum_kg(H1B.ov(hB,pB),T3C.PPpHhH,'me,ABeImJ->ABIJ')...
+einsum_kg(H1B.ov(hB,PB),T3C.PPPHhH,'mE,AEBImJ->ABIJ')...
-einsum_kg(H1B.ov(HB,pB),T3C.PPpHHH,'Me,ABeIMJ->ABIJ')...
+einsum_kg(H1B.ov(HB,PB),T3C.PPPHHH,'ME,AEBIMJ->ABIJ');

    X2B_PPHH = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10;
    X2B(PA,PB,HA,HB) = X2B(PA,PB,HA,HB) + X2B_PPHH;

    
end