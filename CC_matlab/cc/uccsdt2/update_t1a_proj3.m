function [X1A] = update_t1a_proj3(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys,shift)

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


D1 = +0.25*einsum_kg(sys.vA_oovv(hA,hA,pA,pA),T3A.Ppphhh,'mnef,Aefimn->Ai')...
+0.5*einsum_kg(sys.vA_oovv(hA,hA,PA,pA),T3A.PPphhh,'mnEf,AEfimn->Ai')...
+0.25*einsum_kg(sys.vA_oovv(hA,hA,PA,PA),T3A.PPPhhh,'mnEF,AEFimn->Ai')...
-0.5*einsum_kg(sys.vA_oovv(HA,hA,pA,pA),T3A.PpphhH,'Mnef,AefinM->Ai')...
-1*einsum_kg(sys.vA_oovv(HA,hA,PA,pA),T3A.PPphhH,'MnEf,AEfinM->Ai')...
-0.5*einsum_kg(sys.vA_oovv(HA,hA,PA,PA),T3A.PPPhhH,'MnEF,AEFinM->Ai')...
+0.25*einsum_kg(sys.vA_oovv(HA,HA,pA,pA),T3A.PpphHH,'MNef,AefiMN->Ai')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,PA,pA),T3A.PPphHH,'MNEf,AEfiMN->Ai')...
+0.25*einsum_kg(sys.vA_oovv(HA,HA,PA,PA),T3A.PPPhHH,'MNEF,AEFiMN->Ai');

D2 = +1*einsum_kg(sys.vB_oovv(hA,hB,pA,pB),T3B.Ppphhh,'mnef,Aefimn->Ai')...
+1*einsum_kg(sys.vB_oovv(hA,hB,pA,PB),T3B.PpPhhh,'mneF,AeFimn->Ai')...
+1*einsum_kg(sys.vB_oovv(hA,hB,PA,pB),T3B.PPphhh,'mnEf,AEfimn->Ai')...
+1*einsum_kg(sys.vB_oovv(hA,hB,PA,PB),T3B.PPPhhh,'mnEF,AEFimn->Ai')...
+1*einsum_kg(sys.vB_oovv(hA,HB,pA,pB),T3B.PpphhH,'mNef,AefimN->Ai')...
+1*einsum_kg(sys.vB_oovv(hA,HB,pA,PB),T3B.PpPhhH,'mNeF,AeFimN->Ai')...
+1*einsum_kg(sys.vB_oovv(hA,HB,PA,pB),T3B.PPphhH,'mNEf,AEfimN->Ai')...
+1*einsum_kg(sys.vB_oovv(hA,HB,PA,PB),T3B.PPPhhH,'mNEF,AEFimN->Ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,pA,pB),T3B.PpphHh,'Mnef,AefiMn->Ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,pA,PB),T3B.PpPhHh,'MneF,AeFiMn->Ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,PA,pB),T3B.PPphHh,'MnEf,AEfiMn->Ai')...
+1*einsum_kg(sys.vB_oovv(HA,hB,PA,PB),T3B.PPPhHh,'MnEF,AEFiMn->Ai')...
+1*einsum_kg(sys.vB_oovv(HA,HB,pA,pB),T3B.PpphHH,'MNef,AefiMN->Ai')...
+1*einsum_kg(sys.vB_oovv(HA,HB,pA,PB),T3B.PpPhHH,'MNeF,AeFiMN->Ai')...
+1*einsum_kg(sys.vB_oovv(HA,HB,PA,pB),T3B.PPphHH,'MNEf,AEfiMN->Ai')...
+1*einsum_kg(sys.vB_oovv(HA,HB,PA,PB),T3B.PPPhHH,'MNEF,AEFiMN->Ai');

D3 = +0.25*einsum_kg(sys.vC_oovv(hB,hB,pB,pB),T3C.Ppphhh,'mnef,Aefimn->Ai')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,PB,pB),T3C.PPphhh,'mnEf,AEfimn->Ai')...
+0.25*einsum_kg(sys.vC_oovv(hB,hB,PB,PB),T3C.PPPhhh,'mnEF,AEFimn->Ai')...
-0.5*einsum_kg(sys.vC_oovv(HB,hB,pB,pB),T3C.PpphhH,'Mnef,AefinM->Ai')...
-1*einsum_kg(sys.vC_oovv(HB,hB,PB,pB),T3C.PPphhH,'MnEf,AEfinM->Ai')...
-0.5*einsum_kg(sys.vC_oovv(HB,hB,PB,PB),T3C.PPPhhH,'MnEF,AEFinM->Ai')...
+0.25*einsum_kg(sys.vC_oovv(HB,HB,pB,pB),T3C.PpphHH,'MNef,AefiMN->Ai')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,pB),T3C.PPphHH,'MNEf,AEfiMN->Ai')...
+0.25*einsum_kg(sys.vC_oovv(HB,HB,PB,PB),T3C.PPPhHH,'MNEF,AEFiMN->Ai');


    X1A_Ph = D1 + D2 + D3;

    X1A(PA,hA) = X1A(PA,hA) + X1A_Ph;

end