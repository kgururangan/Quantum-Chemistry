function [X2A] = update_t2a_proj2(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys,shift)

    % <ijab|H(CCS)|0>
    I2A_vooo = sys.vA_vooo ...
              -0.5*einsum_kg(sys.vA_oooo,t1a,'nmij,an->amij')...
              +einsum_kg(einsum_kg(sys.vA_vovv,t1a,'amef,ei->amif'),t1a,'amif,fj->amij')...
              +einsum_kg(sys.vA_voov,t1a,'amie,ej->amij')...
                   -einsum_kg(sys.vA_voov,t1a,'amje,ei->amij')...
              -0.5*einsum_kg(einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmef,fj->nmej'),t1a,'nmej,an->amej'),t1a,'amej,ei->amij');

    I2A_vvov  = sys.vA_vvov ...
                +0.5*einsum_kg(sys.vA_vvvv,t1a,'abfe,fi->abie')...
                +einsum_kg(einsum_kg(sys.vA_ooov,t1a,'mnie,am->anie'),t1a,'anie,bn->abie');
    
    % H(CCS) intermediates
    h1A_ov = sys.fa_ov...
             +einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me')...
             +einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'); 
    
    h1A_vv = sys.fa_vv...
             -einsum_kg(h1A_ov,t1a,'me,am->ae')...
             +einsum_kg(sys.vA_vovv,t1a,'amef,fm->ae')...
             +einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae');
     
    h1A_oo = sys.fa_oo...
             +einsum_kg(h1A_ov,t1a,'me,ei->mi')...
             +einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi')...
             +einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi');
     
    h2A_oooo = sys.vA_oooo...
                +einsum_kg(sys.vA_oovo,t1a,'mnej,ei->mnij')...
	     			-einsum_kg(sys.vA_oovo,t1a,'mnei,ej->mnij')...
                +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ei->mnif'),t1a,'mnif,fj->mnij');
            
    h2A_vvvv = sys.vA_vvvv...
                -einsum_kg(sys.vA_ovvv,t1a,'mbef,am->abef')...
	     			+einsum_kg(sys.vA_ovvv,t1a,'maef,bm->abef')...
                +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bn->mbef'),t1a,'mbef,am->abef');
            
    h2A_voov = sys.vA_voov ...
               -einsum_kg(sys.vA_ooov,t1a,'nmie,an->amie')...
               +einsum_kg(sys.vA_vovv,t1a,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie');
            
    h2B_voov = sys.vB_voov ...
               -einsum_kg(sys.vB_ooov,t1a,'nmie,an->amie')...
               +einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie');
    
    % VT2 intermediates
    VT2_1A_oo = 0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi')...
                +einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi');
            
    VT2_1A_vv = -0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,afmn->ae')...
                -einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ae');
            
    % note weights of 0.5 and 1 to account on (V T2A^2) and (V T2A T2B),
    % respectively to account for A(ij)A(ab) in contraction term
    VT2_2A_voov = 0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,afin->amie')...
                  +einsum_kg(sys.vB_oovv,t2b,'mnef,afin->amie');  
    
    % weight of 0.5 on (V T2B^2) term to correct for A(ij)A(ab) in
    % contraction term
    VT2_2B_voov = 0.5*einsum_kg(sys.vC_oovv,t2b,'mnef,afin->amie');
    
    VT2_2A_oooo = 0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efij->mnij');
    
    % contraction intermdiates
    I1A_oo = h1A_oo + VT2_1A_oo;
    I1A_vv = h1A_vv + VT2_1A_vv;
    I2A_voov = h2A_voov + VT2_2A_voov;
    I2B_voov = h2B_voov + VT2_2B_voov;
    I2A_oooo = h2A_oooo + VT2_2A_oooo;
    
    % <ijab| H(CCS) | 0> + <ijab| (H(CCS)T2)_C |0> + <ijab| 0.5*(H(CCS)T2^2)_C |0>
    D1 = -einsum_kg(I2A_vooo,t1a,'amij,bm->abij'); 
    
    D2 = einsum_kg(I2A_vvov,t1a,'abie,ej->abij');
    
    D3 = einsum_kg(I1A_vv,t2a,'ae,ebij->abij');
    
    D4 = -einsum_kg(I1A_oo,t2a,'mi,abmj->abij');
    
    D5 = einsum_kg(I2A_voov,t2a,'amie,ebmj->abij');
    
    D6 = einsum_kg(I2B_voov,t2b,'amie,bejm->abij');
    
    D7 = 0.5*einsum_kg(h2A_vvvv,t2a,'abef,efij->abij');
    
    D8 = 0.5*einsum_kg(I2A_oooo,t2a,'mnij,abmn->abij');
    
    % diagrams that have A(ab)
    D13 = D1 + D3;
    D13 = D13 - permute(D13,[2,1,3,4]);
    
    % diagrams that have A(ij)
    D24 = D2 + D4;
    D24 = D24 - permute(D24,[1,2,4,3]);
        
    % diagrams that have A(ab)A(ij)
    D56 = D5 + D6;
    D56 = D56 - permute(D56,[2,1,3,4]) - permute(D56,[1,2,4,3]) + permute(D56,[2,1,4,3]);
    
    % total contribution
    X2A = sys.vA_vvoo + D13 + D24 + D56 + D7 + D8;
    
    % CCSDT part    
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

    h1B_ov = sys.fb_ov + einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me') + einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
    h2B_ooov = sys.vB_ooov + einsum_kg(sys.vB_oovv,t1a,'mnfe,fi->mnie');
    h2A_ooov = sys.vA_ooov + einsum_kg(sys.vA_oovv,t1a,'mnfe,fi->mnie');
    h2A_vovv = sys.vA_vovv - einsum_kg(sys.vA_oovv,t1a,'mnfe,an->amef');
    h2B_vovv = sys.vB_vovv - einsum_kg(sys.vB_oovv,t1a,'nmef,an->amef');

    H1A.oo = h1A_oo; H1A.vv = h1A_vv; H1A.ov = h1A_ov;
    H1B.ov = h1B_ov; 
    H2A.ooov = h2A_ooov; H2A.vovv = h2A_vovv;
    H2B.ooov = h2B_ooov; H2B.vovv = h2B_vovv;
    



D1 = +einsum_kg(H1A.ov(hA,pA),T3A.PpphHH,'me,AbemIJ->AbIJ')...
-einsum_kg(H1A.ov(hA,PA),T3A.PPphHH,'mE,AEbmIJ->AbIJ')...
+einsum_kg(H1A.ov(HA,pA),T3A.PppHHH,'Me,AbeIJM->AbIJ')...
-einsum_kg(H1A.ov(HA,PA),T3A.PPpHHH,'ME,AEbIJM->AbIJ');

D2 = +einsum_kg(H1B.ov(hB,pB),T3B.PppHHh,'me,AbeIJm->AbIJ')...
+einsum_kg(H1B.ov(hB,PB),T3B.PpPHHh,'mE,AbEIJm->AbIJ')...
+einsum_kg(H1B.ov(HB,pB),T3B.PppHHH,'Me,AbeIJM->AbIJ')...
+einsum_kg(H1B.ov(HB,PB),T3B.PpPHHH,'ME,AbEIJM->AbIJ');

D3 = -einsum_kg(H2B.ooov(hA,hB,HA,pB),T3B.PpphHh,'mnIf,AbfmJn->AbIJ')...
-einsum_kg(H2B.ooov(hA,hB,HA,PB),T3B.PpPhHh,'mnIF,AbFmJn->AbIJ')...
-einsum_kg(H2B.ooov(hA,HB,HA,pB),T3B.PpphHH,'mNIf,AbfmJN->AbIJ')...
-einsum_kg(H2B.ooov(hA,HB,HA,PB),T3B.PpPhHH,'mNIF,AbFmJN->AbIJ')...
-einsum_kg(H2B.ooov(HA,hB,HA,pB),T3B.PppHHh,'MnIf,AbfMJn->AbIJ')...
-einsum_kg(H2B.ooov(HA,hB,HA,PB),T3B.PpPHHh,'MnIF,AbFMJn->AbIJ')...
-einsum_kg(H2B.ooov(HA,HB,HA,pB),T3B.PppHHH,'MNIf,AbfMJN->AbIJ')...
-einsum_kg(H2B.ooov(HA,HB,HA,PB),T3B.PpPHHH,'MNIF,AbFMJN->AbIJ');
D3 = D3 - permute(D3,[1,2,4,3]);

D4 = +0.5*einsum_kg(H2A.ooov(hA,hA,HA,pA),T3A.PpphhH,'mnIf,AbfmnJ->AbIJ')...
-0.5*einsum_kg(H2A.ooov(hA,hA,HA,PA),T3A.PPphhH,'mnIF,AFbmnJ->AbIJ')...
-einsum_kg(H2A.ooov(HA,hA,HA,pA),T3A.PpphHH,'MnIf,AbfnMJ->AbIJ')...
+einsum_kg(H2A.ooov(HA,hA,HA,PA),T3A.PPphHH,'MnIF,AFbnMJ->AbIJ')...
-0.5*einsum_kg(H2A.ooov(HA,HA,HA,pA),T3A.PppHHH,'MNIf,AbfMJN->AbIJ')...
+0.5*einsum_kg(H2A.ooov(HA,HA,HA,PA),T3A.PPpHHH,'MNIF,AFbMJN->AbIJ');
D4 = D4 - permute(D4,[1,2,4,3]);

D5 = +0.5*einsum_kg(H2A.vovv(PA,hA,pA,pA),T3A.ppphHH,'Anef,ebfnIJ->AbIJ')...
+einsum_kg(H2A.vovv(PA,hA,PA,pA),T3A.PpphHH,'AnEf,EbfnIJ->AbIJ')...
-0.5*einsum_kg(H2A.vovv(PA,hA,PA,PA),T3A.PPphHH,'AnEF,EFbnIJ->AbIJ')...
+0.5*einsum_kg(H2A.vovv(PA,HA,pA,pA),T3A.pppHHH,'ANef,ebfIJN->AbIJ')...
+einsum_kg(H2A.vovv(PA,HA,PA,pA),T3A.PppHHH,'ANEf,EbfIJN->AbIJ')...
-0.5*einsum_kg(H2A.vovv(PA,HA,PA,PA),T3A.PPpHHH,'ANEF,EFbIJN->AbIJ')...
+0.5*einsum_kg(H2A.vovv(pA,hA,pA,pA),T3A.PpphHH,'bnef,AefnIJ->AbIJ')...
-einsum_kg(H2A.vovv(pA,hA,PA,pA),T3A.PPphHH,'bnEf,EAfnIJ->AbIJ')...
-0.5*einsum_kg(H2A.vovv(pA,hA,PA,PA),T3A.PPPhHH,'bnEF,EAFnIJ->AbIJ')...
+0.5*einsum_kg(H2A.vovv(pA,HA,pA,pA),T3A.PppHHH,'bNef,AefIJN->AbIJ')...
-einsum_kg(H2A.vovv(pA,HA,PA,pA),T3A.PPpHHH,'bNEf,EAfIJN->AbIJ')...
-0.5*einsum_kg(H2A.vovv(pA,HA,PA,PA),T3A.PPPHHH,'bNEF,EAFIJN->AbIJ');

D6 = +einsum_kg(H2B.vovv(PA,hB,pA,pB),T3B.pppHHh,'Anef,ebfIJn->AbIJ')...
+einsum_kg(H2B.vovv(PA,hB,pA,PB),T3B.ppPHHh,'AneF,ebFIJn->AbIJ')...
+einsum_kg(H2B.vovv(PA,hB,PA,pB),T3B.PppHHh,'AnEf,EbfIJn->AbIJ')...
+einsum_kg(H2B.vovv(PA,hB,PA,PB),T3B.PpPHHh,'AnEF,EbFIJn->AbIJ')...
+einsum_kg(H2B.vovv(PA,HB,pA,pB),T3B.pppHHH,'ANef,ebfIJN->AbIJ')...
+einsum_kg(H2B.vovv(PA,HB,pA,PB),T3B.ppPHHH,'ANeF,ebFIJN->AbIJ')...
+einsum_kg(H2B.vovv(PA,HB,PA,pB),T3B.PppHHH,'ANEf,EbfIJN->AbIJ')...
+einsum_kg(H2B.vovv(PA,HB,PA,PB),T3B.PpPHHH,'ANEF,EbFIJN->AbIJ')...
+einsum_kg(H2B.vovv(pA,hB,pA,pB),T3B.PppHHh,'bnef,AefIJn->AbIJ')...
+einsum_kg(H2B.vovv(pA,hB,pA,PB),T3B.PpPHHh,'bneF,AeFIJn->AbIJ')...
-einsum_kg(H2B.vovv(pA,hB,PA,pB),T3B.PPpHHh,'bnEf,EAfIJn->AbIJ')...
-einsum_kg(H2B.vovv(pA,hB,PA,PB),T3B.PPPHHh,'bnEF,EAFIJn->AbIJ')...
+einsum_kg(H2B.vovv(pA,HB,pA,pB),T3B.PppHHH,'bNef,AefIJN->AbIJ')...
+einsum_kg(H2B.vovv(pA,HB,pA,PB),T3B.PpPHHH,'bNeF,AeFIJN->AbIJ')...
-einsum_kg(H2B.vovv(pA,HB,PA,pB),T3B.PPpHHH,'bNEf,EAfIJN->AbIJ')...
-einsum_kg(H2B.vovv(pA,HB,PA,PB),T3B.PPPHHH,'bNEF,EAFIJN->AbIJ');
    
    X2A_PpHH = ( D1 + D2 + D3 + D4 + D5 + D6 );

    X2A(PA,pA,HA,HA) = X2A(PA,pA,HA,HA) + X2A_PpHH;


end