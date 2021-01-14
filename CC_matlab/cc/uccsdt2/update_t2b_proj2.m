function [X2B] = update_t2b_proj2(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys,shift)

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
    


D1 = -0.5*einsum_kg(H2A.ooov(hA,hA,HA,pA),T3B.PpphhH,'mnIf,AfbmnJ->AbIJ')...
-0.5*einsum_kg(H2A.ooov(hA,hA,HA,PA),T3B.PPphhH,'mnIF,AFbmnJ->AbIJ')...
+einsum_kg(H2A.ooov(HA,hA,HA,pA),T3B.PpphHH,'MnIf,AfbnMJ->AbIJ')...
+einsum_kg(H2A.ooov(HA,hA,HA,PA),T3B.PPphHH,'MnIF,AFbnMJ->AbIJ')...
-0.5*einsum_kg(H2A.ooov(HA,HA,HA,pA),T3B.PppHHH,'MNIf,AfbMNJ->AbIJ')...
-0.5*einsum_kg(H2A.ooov(HA,HA,HA,PA),T3B.PPpHHH,'MNIF,AFbMNJ->AbIJ');

D2 = +einsum_kg(H2B.oovo(hA,hB,pA,HB),T3B.PpphHh,'nmfJ,AfbnIm->AbIJ')...
+einsum_kg(H2B.oovo(hA,hB,PA,HB),T3B.PPphHh,'nmFJ,AFbnIm->AbIJ')...
+einsum_kg(H2B.oovo(hA,HB,pA,HB),T3B.PpphHH,'nMfJ,AfbnIM->AbIJ')...
+einsum_kg(H2B.oovo(hA,HB,PA,HB),T3B.PPphHH,'nMFJ,AFbnIM->AbIJ')...
-einsum_kg(H2B.oovo(HA,hB,pA,HB),T3B.PppHHh,'NmfJ,AfbINm->AbIJ')...
-einsum_kg(H2B.oovo(HA,hB,PA,HB),T3B.PPpHHh,'NmFJ,AFbINm->AbIJ')...
-einsum_kg(H2B.oovo(HA,HB,pA,HB),T3B.PppHHH,'NMfJ,AfbINM->AbIJ')...
-einsum_kg(H2B.oovo(HA,HB,PA,HB),T3B.PPpHHH,'NMFJ,AFbINM->AbIJ');

D3 = -0.5*einsum_kg(H2C.ooov(hB,hB,HB,pB),T3C.PppHhh,'mnJf,AfbInm->AbIJ')...
-0.5*einsum_kg(H2C.ooov(hB,hB,HB,PB),T3C.PPpHhh,'mnJF,AFbInm->AbIJ')...
-einsum_kg(H2C.ooov(HB,hB,HB,pB),T3C.PppHhH,'MnJf,AfbInM->AbIJ')...
-einsum_kg(H2C.ooov(HB,hB,HB,PB),T3C.PPpHhH,'MnJF,AFbInM->AbIJ')...
-0.5*einsum_kg(H2C.ooov(HB,HB,HB,pB),T3C.PppHHH,'MNJf,AfbINM->AbIJ')...
-0.5*einsum_kg(H2C.ooov(HB,HB,HB,PB),T3C.PPpHHH,'MNJF,AFbINM->AbIJ');

D4 = -einsum_kg(H2B.ooov(hA,hB,HA,pB),T3C.PpphhH,'mnIf,AfbmnJ->AbIJ')...
-einsum_kg(H2B.ooov(hA,hB,HA,PB),T3C.PPphhH,'mnIF,AFbmnJ->AbIJ')...
-einsum_kg(H2B.ooov(hA,HB,HA,pB),T3C.PpphHH,'mNIf,AfbmNJ->AbIJ')...
-einsum_kg(H2B.ooov(hA,HB,HA,PB),T3C.PPphHH,'mNIF,AFbmNJ->AbIJ')...
-einsum_kg(H2B.ooov(HA,hB,HA,pB),T3C.PppHhH,'MnIf,AfbMnJ->AbIJ')...
-einsum_kg(H2B.ooov(HA,hB,HA,PB),T3C.PPpHhH,'MnIF,AFbMnJ->AbIJ')...
-einsum_kg(H2B.ooov(HA,HB,HA,pB),T3C.PppHHH,'MNIf,AfbMNJ->AbIJ')...
-einsum_kg(H2B.ooov(HA,HB,HA,PB),T3C.PPpHHH,'MNIF,AFbMNJ->AbIJ');

D5 = -0.5*einsum_kg(H2A.vovv(PA,hA,pA,pA),T3B.ppphHH,'Anef,efbnIJ->AbIJ')...
-einsum_kg(H2A.vovv(PA,hA,PA,pA),T3B.PpphHH,'AnEf,EfbnIJ->AbIJ')...
-0.5*einsum_kg(H2A.vovv(PA,hA,PA,PA),T3B.PPphHH,'AnEF,EFbnIJ->AbIJ')...
+0.5*einsum_kg(H2A.vovv(PA,HA,pA,pA),T3B.pppHHH,'ANef,efbINJ->AbIJ')...
+einsum_kg(H2A.vovv(PA,HA,PA,pA),T3B.PppHHH,'ANEf,EfbINJ->AbIJ')...
+0.5*einsum_kg(H2A.vovv(PA,HA,PA,PA),T3B.PPpHHH,'ANEF,EFbINJ->AbIJ');

D6 = +einsum_kg(H2B.vovv(PA,hB,pA,pB),T3C.pppHhH,'Anef,efbInJ->AbIJ')...
+einsum_kg(H2B.vovv(PA,hB,pA,PB),T3C.pPpHhH,'AneF,eFbInJ->AbIJ')...
+einsum_kg(H2B.vovv(PA,hB,PA,pB),T3C.PppHhH,'AnEf,EfbInJ->AbIJ')...
+einsum_kg(H2B.vovv(PA,hB,PA,PB),T3C.PPpHhH,'AnEF,EFbInJ->AbIJ')...
+einsum_kg(H2B.vovv(PA,HB,pA,pB),T3C.pppHHH,'ANef,efbINJ->AbIJ')...
+einsum_kg(H2B.vovv(PA,HB,pA,PB),T3C.pPpHHH,'ANeF,eFbINJ->AbIJ')...
+einsum_kg(H2B.vovv(PA,HB,PA,pB),T3C.PppHHH,'ANEf,EfbINJ->AbIJ')...
+einsum_kg(H2B.vovv(PA,HB,PA,PB),T3C.PPpHHH,'ANEF,EFbINJ->AbIJ');

D7 = -einsum_kg(H2B.ovvv(hA,pB,pA,pB),T3B.PpphHH,'nbfe,AfenIJ->AbIJ')...
-einsum_kg(H2B.ovvv(hA,pB,pA,PB),T3B.PpPhHH,'nbfE,AfEnIJ->AbIJ')...
-einsum_kg(H2B.ovvv(hA,pB,PA,pB),T3B.PPphHH,'nbFe,AFenIJ->AbIJ')...
-einsum_kg(H2B.ovvv(hA,pB,PA,PB),T3B.PPPhHH,'nbFE,AFEnIJ->AbIJ')...
+einsum_kg(H2B.ovvv(HA,pB,pA,pB),T3B.PppHHH,'Nbfe,AfeINJ->AbIJ')...
+einsum_kg(H2B.ovvv(HA,pB,pA,PB),T3B.PpPHHH,'NbfE,AfEINJ->AbIJ')...
+einsum_kg(H2B.ovvv(HA,pB,PA,pB),T3B.PPpHHH,'NbFe,AFeINJ->AbIJ')...
+einsum_kg(H2B.ovvv(HA,pB,PA,PB),T3B.PPPHHH,'NbFE,AFEINJ->AbIJ');

D8 = +0.5*einsum_kg(H2C.vovv(pB,hB,pB,pB),T3C.PppHhH,'bnef,AfeInJ->AbIJ')...
-einsum_kg(H2C.vovv(pB,hB,PB,pB),T3C.PPpHhH,'bnEf,AEfInJ->AbIJ')...
+0.5*einsum_kg(H2C.vovv(pB,hB,PB,PB),T3C.PPPHhH,'bnEF,AFEInJ->AbIJ')...
+0.5*einsum_kg(H2C.vovv(pB,HB,pB,pB),T3C.PppHHH,'bNef,AfeINJ->AbIJ')...
-einsum_kg(H2C.vovv(pB,HB,PB,pB),T3C.PPpHHH,'bNEf,AEfINJ->AbIJ')...
+0.5*einsum_kg(H2C.vovv(pB,HB,PB,PB),T3C.PPPHHH,'bNEF,AFEINJ->AbIJ');

D9 = -einsum_kg(H1A.ov(hA,pA),T3B.PpphHH,'me,AebmIJ->AbIJ')...
-einsum_kg(H1A.ov(hA,PA),T3B.PPphHH,'mE,AEbmIJ->AbIJ')...
+einsum_kg(H1A.ov(HA,pA),T3B.PppHHH,'Me,AebIMJ->AbIJ')...
+einsum_kg(H1A.ov(HA,PA),T3B.PPpHHH,'ME,AEbIMJ->AbIJ');

D10 = +einsum_kg(H1B.ov(hB,pB),T3C.PppHhH,'me,AebImJ->AbIJ')...
+einsum_kg(H1B.ov(hB,PB),T3C.PPpHhH,'mE,AEbImJ->AbIJ')...
+einsum_kg(H1B.ov(HB,pB),T3C.PppHHH,'Me,AebIMJ->AbIJ')...
+einsum_kg(H1B.ov(HB,PB),T3C.PPpHHH,'ME,AEbIMJ->AbIJ');

    X2B_PpHH = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10;
    X2B(PA,pB,HA,HB) = X2B(PA,pB,HA,HB) + X2B_PpHH;

    
end