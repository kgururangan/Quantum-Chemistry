function [X1A] = update_t1a_proj4(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys,shift)

    chi1A_vv = sys.fa_vv + einsum_kg(sys.vA_vovv,t1a,'anef,fn->ae') + einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae');
    chi1A_oo = sys.fa_oo + einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi') + einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi');
    h1A_ov = sys.fa_ov + einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me') + einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');
    h1B_ov = sys.fb_ov + einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me') + einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
    h1A_oo = chi1A_oo + einsum_kg(h1A_ov,t1a,'me,ei->mi');
       
    M11 = sys.fa_vo - einsum_kg(h1A_oo,t1a,'mi,am->ai') + einsum_kg(chi1A_vv,t1a,'ae,ei->ai')...
          +einsum_kg(sys.vA_voov,t1a,'anif,fn->ai') + einsum_kg(sys.vB_voov,t1b,'anif,fn->ai');
      
    h2A_ooov = sys.vA_ooov + einsum_kg(sys.vA_oovv,t1a,'mnfe,fi->mnie'); 
    h2B_ooov = sys.vB_ooov + einsum_kg(sys.vB_oovv,t1a,'mnfe,fi->mnie');
    h2A_vovv = sys.vA_vovv - einsum_kg(sys.vA_oovv,t1a,'mnfe,an->amef');
    h2B_vovv = sys.vB_vovv - einsum_kg(sys.vB_oovv,t1a,'nmef,an->amef');
      
    CCS_T2 = +einsum_kg(h1A_ov,t2a,'me,aeim->ai') + einsum_kg(h1B_ov,t2b,'me,aeim->ai')...
             -0.5*einsum_kg(h2A_ooov,t2a,'mnif,afmn->ai') - einsum_kg(h2B_ooov,t2b,'mnif,afmn->ai')...
             +0.5*einsum_kg(h2A_vovv,t2a,'anef,efin->ai') + einsum_kg(h2B_vovv,t2b,'anef,efin->ai');
       
    X1A = M11 + CCS_T2; 
    
    % CCSDT part
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;


D1 = +0.25*einsum_kg(sys.vA_oovv(hA,hA,pA,pA),T3A.ppphhh,'mnef,aefimn->ai')...
-0.5*einsum_kg(sys.vA_oovv(hA,hA,PA,pA),T3A.Ppphhh,'mnEf,Eafimn->ai')...
+0.25*einsum_kg(sys.vA_oovv(hA,hA,PA,PA),T3A.PPphhh,'mnEF,EFaimn->ai')...
-0.5*einsum_kg(sys.vA_oovv(HA,hA,pA,pA),T3A.ppphhH,'Mnef,aefinM->ai')...
+1*einsum_kg(sys.vA_oovv(HA,hA,PA,pA),T3A.PpphhH,'MnEf,EafinM->ai')...
-0.5*einsum_kg(sys.vA_oovv(HA,hA,PA,PA),T3A.PPphhH,'MnEF,EFainM->ai')...
+0.25*einsum_kg(sys.vA_oovv(HA,HA,pA,pA),T3A.ppphHH,'MNef,aefiMN->ai')...
-0.5*einsum_kg(sys.vA_oovv(HA,HA,PA,pA),T3A.PpphHH,'MNEf,EafiMN->ai')...
+0.25*einsum_kg(sys.vA_oovv(HA,HA,PA,PA),T3A.PPphHH,'MNEF,EFaiMN->ai');

D2 = +1*einsum_kg(sys.vB_oovv(hA,hB,pA,pB),T3B.ppphhh,'mnef,aefimn->ai')...
+1*einsum_kg(sys.vB_oovv(hA,hB,pA,PB),T3B.ppPhhh,'mneF,aeFimn->ai')...
-1*einsum_kg(sys.vB_oovv(hA,hB,PA,pB),T3B.Ppphhh,'mnEf,Eafimn->ai')...
-1*einsum_kg(sys.vB_oovv(hA,hB,PA,PB),T3B.PpPhhh,'mnEF,EaFimn->ai')...
+1*einsum_kg(sys.vB_oovv(hA,HB,pA,pB),T3B.ppphhH,'mNef,aefimN->ai')...
+1*einsum_kg(sys.vB_oovv(hA,HB,pA,PB),T3B.ppPhhH,'mNeF,aeFimN->ai')...
-1*einsum_kg(sys.vB_oovv(hA,HB,PA,pB),T3B.PpphhH,'mNEf,EafimN->ai')...
-1*einsum_kg(sys.vB_oovv(hA,HB,PA,PB),T3B.PpPhhH,'mNEF,EaFimN->ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,pA,pB),T3B.ppphHh,'Mnef,aefiMn->ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,pA,PB),T3B.ppPhHh,'MneF,aeFiMn->ai')...
-1*einsum_kg(sys.vB_oovv(HA,hB,PA,pB),T3B.PpphHh,'MnEf,EafiMn->ai')...
-1*einsum_kg(sys.vB_oovv(HA,hB,PA,PB),T3B.PpPhHh,'MnEF,EaFiMn->ai')...
+1*einsum_kg(sys.vB_oovv(HA,HB,pA,pB),T3B.ppphHH,'MNef,aefiMN->ai')...
+1*einsum_kg(sys.vB_oovv(HA,HB,pA,PB),T3B.ppPhHH,'MNeF,aeFiMN->ai')...
-1*einsum_kg(sys.vB_oovv(HA,HB,PA,pB),T3B.PpphHH,'MNEf,EafiMN->ai')...
-1*einsum_kg(sys.vB_oovv(HA,HB,PA,PB),T3B.PpPhHH,'MNEF,EaFiMN->ai');

D3 = +0.25*einsum_kg(sys.vC_oovv(hB,hB,pB,pB),T3C.ppphhh,'mnef,aefimn->ai')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,PB,pB),T3C.pPphhh,'mnEf,aEfimn->ai')...
+0.25*einsum_kg(sys.vC_oovv(hB,hB,PB,PB),T3C.pPPhhh,'mnEF,aEFimn->ai')...
-0.5*einsum_kg(sys.vC_oovv(HB,hB,pB,pB),T3C.ppphhH,'Mnef,aefinM->ai')...
-1*einsum_kg(sys.vC_oovv(HB,hB,PB,pB),T3C.pPphhH,'MnEf,aEfinM->ai')...
-0.5*einsum_kg(sys.vC_oovv(HB,hB,PB,PB),T3C.pPPhhH,'MnEF,aEFinM->ai')...
+0.25*einsum_kg(sys.vC_oovv(HB,HB,pB,pB),T3C.ppphHH,'MNef,aefiMN->ai')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,pB),T3C.pPphHH,'MNEf,aEfiMN->ai')...
+0.25*einsum_kg(sys.vC_oovv(HB,HB,PB,PB),T3C.pPPhHH,'MNEF,aEFiMN->ai');


    X1A_ph = D1 + D2 + D3;

    X1A(pA,hA) = X1A(pA,hA) + X1A_ph;

end