function [X1B] = update_t1b_proj4(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys,shift)

    chi1B_vv = sys.fb_vv + einsum_kg(sys.vC_vovv,t1b,'anef,fn->ae') + einsum_kg(sys.vB_ovvv,t1a,'nafe,fn->ae');
    chi1B_oo = sys.fb_oo + einsum_kg(sys.vC_ooov,t1b,'mnif,fn->mi') + einsum_kg(sys.vB_oovo,t1b,'nmfi,fn->mi');
    h1A_ov = sys.fa_ov + einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me') + einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');
    h1B_ov = sys.fb_ov + einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me') + einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
    h1B_oo = chi1B_oo + einsum_kg(h1B_ov,t1b,'me,ei->mi');
       
    M11 = sys.fb_vo - einsum_kg(h1B_oo,t1b,'mi,am->ai') + einsum_kg(chi1B_vv,t1b,'ae,ei->ai')...
          +einsum_kg(sys.vC_voov,t1b,'anif,fn->ai') + einsum_kg(sys.vB_ovvo,t1a,'nafi,fn->ai');
      
    h2C_ooov = sys.vC_ooov + einsum_kg(sys.vC_oovv,t1b,'mnfe,fi->mnie'); 
    h2B_oovo = sys.vB_oovo + einsum_kg(sys.vB_oovv,t1b,'nmef,fi->nmei');
    h2C_vovv = sys.vC_vovv - einsum_kg(sys.vC_oovv,t1b,'mnfe,an->amef');
    h2B_ovvv = sys.vB_ovvv - einsum_kg(sys.vB_oovv,t1b,'mnfe,an->mafe');
      
    CCS_T2 = +einsum_kg(h1A_ov,t2b,'me,eami->ai') + einsum_kg(h1B_ov,t2c,'me,aeim->ai')...
             -0.5*einsum_kg(h2C_ooov,t2c,'mnif,afmn->ai') - einsum_kg(h2B_oovo,t2b,'nmfi,fanm->ai')...
             +0.5*einsum_kg(h2C_vovv,t2c,'anef,efin->ai') + einsum_kg(h2B_ovvv,t2b,'nafe,feni->ai');
       
    X1B = M11 + CCS_T2; 
    
    % CCSDT part
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;


D1 = +0.25*einsum_kg(sys.vC_oovv(hB,hB,pB,pB),T3D.ppphhh,'mnef,aefimn->ai')...
-0.5*einsum_kg(sys.vC_oovv(hB,hB,PB,pB),T3D.Ppphhh,'mnEf,Eafimn->ai')...
+0.25*einsum_kg(sys.vC_oovv(hB,hB,PB,PB),T3D.PPphhh,'mnEF,EFaimn->ai')...
-0.5*einsum_kg(sys.vC_oovv(HB,hB,pB,pB),T3D.ppphhH,'Mnef,aefinM->ai')...
+1*einsum_kg(sys.vC_oovv(HB,hB,PB,pB),T3D.PpphhH,'MnEf,EafinM->ai')...
-0.5*einsum_kg(sys.vC_oovv(HB,hB,PB,PB),T3D.PPphhH,'MnEF,EFainM->ai')...
+0.25*einsum_kg(sys.vC_oovv(HB,HB,pB,pB),T3D.ppphHH,'MNef,aefiMN->ai')...
-0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,pB),T3D.PpphHH,'MNEf,EafiMN->ai')...
+0.25*einsum_kg(sys.vC_oovv(HB,HB,PB,PB),T3D.PPphHH,'MNEF,EFaiMN->ai');

D2 = +1*einsum_kg(sys.vB_oovv(hA,hB,pA,pB),T3C.ppphhh,'mnef,efamni->ai')...
+1*einsum_kg(sys.vB_oovv(hA,hB,pA,PB),T3C.pPphhh,'mneF,eFamni->ai')...
+1*einsum_kg(sys.vB_oovv(hA,hB,PA,pB),T3C.Ppphhh,'mnEf,Efamni->ai')...
+1*einsum_kg(sys.vB_oovv(hA,hB,PA,PB),T3C.PPphhh,'mnEF,EFamni->ai')...
-1*einsum_kg(sys.vB_oovv(hA,HB,pA,pB),T3C.ppphhH,'mNef,efamiN->ai')...
-1*einsum_kg(sys.vB_oovv(hA,HB,pA,PB),T3C.pPphhH,'mNeF,eFamiN->ai')...
-1*einsum_kg(sys.vB_oovv(hA,HB,PA,pB),T3C.PpphhH,'mNEf,EfamiN->ai')...
-1*einsum_kg(sys.vB_oovv(hA,HB,PA,PB),T3C.PPphhH,'mNEF,EFamiN->ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,pA,pB),T3C.pppHhh,'Mnef,efaMni->ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,pA,PB),T3C.pPpHhh,'MneF,eFaMni->ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,PA,pB),T3C.PppHhh,'MnEf,EfaMni->ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,PA,PB),T3C.PPpHhh,'MnEF,EFaMni->ai')...
-1*einsum_kg(sys.vB_oovv(HA,HB,pA,pB),T3C.pppHhH,'MNef,efaMiN->ai')...
-1*einsum_kg(sys.vB_oovv(HA,HB,pA,PB),T3C.pPpHhH,'MNeF,eFaMiN->ai')...
-1*einsum_kg(sys.vB_oovv(HA,HB,PA,pB),T3C.PppHhH,'MNEf,EfaMiN->ai')...
-1*einsum_kg(sys.vB_oovv(HA,HB,PA,PB),T3C.PPpHhH,'MNEF,EFaMiN->ai');

D3 = +0.25*einsum_kg(sys.vA_oovv(hA,hA,pA,pA),T3B.ppphhh,'mnef,efamni->ai')...
+0.5*einsum_kg(sys.vA_oovv(hA,hA,PA,pA),T3B.Ppphhh,'mnEf,Efamni->ai')...
+0.25*einsum_kg(sys.vA_oovv(hA,hA,PA,PA),T3B.PPphhh,'mnEF,EFamni->ai')...
-0.5*einsum_kg(sys.vA_oovv(HA,hA,pA,pA),T3B.ppphHh,'Mnef,efanMi->ai')...
-1*einsum_kg(sys.vA_oovv(HA,hA,PA,pA),T3B.PpphHh,'MnEf,EfanMi->ai')...
-0.5*einsum_kg(sys.vA_oovv(HA,hA,PA,PA),T3B.PPphHh,'MnEF,EFanMi->ai')...
+0.25*einsum_kg(sys.vA_oovv(HA,HA,pA,pA),T3B.pppHHh,'MNef,efaMNi->ai')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,PA,pA),T3B.PppHHh,'MNEf,EfaMNi->ai')...
+0.25*einsum_kg(sys.vA_oovv(HA,HA,PA,PA),T3B.PPpHHh,'MNEF,EFaMNi->ai');
    
    X1B_ph = D1 + D2 + D3;
    X1B(pB,hB) = X1B(pB,hB) + X1B_ph;

end