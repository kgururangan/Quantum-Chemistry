function [X3C_pPPHHH] = update_t3c_proj3(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, sys, shift)
    
    [H1A,H1B,H2A,H2B,H2C,Vt3B,Vt3C] = get_t3c_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys);
    
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
    
    % MM23C
    I2B_vooo = H2B.vooo - einsum_kg(H1B.ov,t2b,'me,aeij->amij');
    I2C_ovoo = H2C.ovoo - einsum_kg(H1B.ov,t2c,'me,ecjk->mcjk');
    I2B_ovoo = H2B.ovoo - einsum_kg(H1A.ov,t2b,'me,ebij->mbij');
    

    M23_D1 = +einsum_kg(H2B.vvov(pA,PB,HA,:) + Vt3B.pPHv,t2c(:,PB,HB,HB),'aBIe,eCJK->aBCIJK');
    
    M23_D2 = -einsum_kg(I2B_vooo(pA,:,HA,HB) + Vt3B.poHH,t2c(PB,PB,:,HB),'amIJ,BCmK->aBCIJK');
  
    M23_D3 = +einsum_kg(H2C.vvvo(PB,PB,:,HB) + Vt3C.PPvH,t2b(pA,:,HA,HB),'BCeK,aeIJ->aBCIJK');

    M23_D4 = -einsum_kg(I2C_ovoo(:,PB,HB,HB) + Vt3C.oPHH,t2b(pA,PB,HA,:),'mCJK,aBIm->aBCIJK');

    M23_D5 = +einsum_kg(H2B.vvvo(pA,PB,:,HB) + Vt3B.pPvH,t2b(:,PB,HA,HB),'aBeJ,eCIK->aBCIJK');

    M23_D6 = -einsum_kg(I2B_ovoo(:,PB,HA,HB) + Vt3B.oPHH,t2b(pA,PB,:,HB),'mBIJ,aCmK->aBCIJK');

    M23_14 = M23_D1 + M23_D4;
    M23_14 = M23_14 - permute(M23_14,[1,3,2,4,5,6]);
    M23_23 = M23_D2 + M23_D3;
    M23_23 = M23_23 - permute(M23_23,[1,2,3,4,6,5]);
    M23_56 = M23_D5 + M23_D6;
    M23_56 = M23_56 - permute(M23_56,[1,3,2,4,5,6]) - permute(M23_56,[1,2,3,4,6,5]) + permute(M23_56,[1,3,2,4,6,5]);
    MM23C = M23_14 + M23_23 + M23_56;

    % (HBar*T3)_C

D1 = -einsum_kg(H1A.oo(hA,HA),T3C.pPPhHH,'mI,aBCmJK->aBCIJK')...
-einsum_kg(H1A.oo(HA,HA),T3C.pPPHHH,'MI,aBCMJK->aBCIJK');

D2 = -einsum_kg(H1B.oo(hB,HB),T3C.pPPHhH,'mJ,aBCImK->aBCIJK')...
-einsum_kg(H1B.oo(HB,HB),T3C.pPPHHH,'MJ,aBCIMK->aBCIJK');
D2 = D2 - permute(D2,[1,2,3,4,6,5]);

D3 = +einsum_kg(H1A.vv(pA,pA),T3C.pPPHHH,'ae,eBCIJK->aBCIJK')...
+einsum_kg(H1A.vv(pA,PA),T3C.PPPHHH,'aE,EBCIJK->aBCIJK');

D4 = -einsum_kg(H1B.vv(PB,pB),T3C.pPpHHH,'Be,aCeIJK->aBCIJK')...
+einsum_kg(H1B.vv(PB,PB),T3C.pPPHHH,'BE,aECIJK->aBCIJK');
D4 = D4 - permute(D4,[1,3,2,4,5,6]);

D5 = +0.5*einsum_kg(H2C.oooo(hB,hB,HB,HB),T3C.pPPHhh,'mnJK,aBCImn->aBCIJK')...
-einsum_kg(H2C.oooo(HB,hB,HB,HB),T3C.pPPHhH,'MnJK,aBCInM->aBCIJK')...
+0.5*einsum_kg(H2C.oooo(HB,HB,HB,HB),T3C.pPPHHH,'MNJK,aBCIMN->aBCIJK');

D6 = +einsum_kg(H2B.oooo(hA,hB,HA,HB),T3C.pPPhhH,'mnIJ,aBCmnK->aBCIJK')...
+einsum_kg(H2B.oooo(hA,HB,HA,HB),T3C.pPPhHH,'mNIJ,aBCmNK->aBCIJK')...
+einsum_kg(H2B.oooo(HA,hB,HA,HB),T3C.pPPHhH,'MnIJ,aBCMnK->aBCIJK')...
+einsum_kg(H2B.oooo(HA,HB,HA,HB),T3C.pPPHHH,'MNIJ,aBCMNK->aBCIJK');
D6 = D6 - permute(D6,[1,2,3,4,6,5]);

D7 = +0.5*einsum_kg(H2C.vvvv(PB,PB,pB,pB),T3C.pppHHH,'BCef,aefIJK->aBCIJK')...
+einsum_kg(H2C.vvvv(PB,PB,PB,pB),T3C.pPpHHH,'BCEf,aEfIJK->aBCIJK')...
+0.5*einsum_kg(H2C.vvvv(PB,PB,PB,PB),T3C.pPPHHH,'BCEF,aEFIJK->aBCIJK');

D8 = -einsum_kg(H2B.vvvv(pA,PB,pA,pB),T3C.pPpHHH,'aBef,eCfIJK->aBCIJK')...
+einsum_kg(H2B.vvvv(pA,PB,pA,PB),T3C.pPPHHH,'aBeF,eFCIJK->aBCIJK')...
-einsum_kg(H2B.vvvv(pA,PB,PA,pB),T3C.PPpHHH,'aBEf,ECfIJK->aBCIJK')...
+einsum_kg(H2B.vvvv(pA,PB,PA,PB),T3C.PPPHHH,'aBEF,EFCIJK->aBCIJK');
D8 = D8 - permute(D8,[1,3,2,4,5,6]);

D9 = +einsum_kg(H2A.voov(pA,hA,HA,pA),T3C.pPPhHH,'amIe,eBCmJK->aBCIJK')...
+einsum_kg(H2A.voov(pA,hA,HA,PA),T3C.PPPhHH,'amIE,EBCmJK->aBCIJK')...
+einsum_kg(H2A.voov(pA,HA,HA,pA),T3C.pPPHHH,'aMIe,eBCMJK->aBCIJK')...
+einsum_kg(H2A.voov(pA,HA,HA,PA),T3C.PPPHHH,'aMIE,EBCMJK->aBCIJK');

D10 = +einsum_kg(H2B.voov(pA,hB,HA,pB),T3D.PPphHH,'amIe,BCemJK->aBCIJK')...
+einsum_kg(H2B.voov(pA,hB,HA,PB),T3D.PPPhHH,'amIE,EBCmJK->aBCIJK')...
+einsum_kg(H2B.voov(pA,HB,HA,pB),T3D.PPpHHH,'aMIe,BCeMJK->aBCIJK')...
+einsum_kg(H2B.voov(pA,HB,HA,PB),T3D.PPPHHH,'aMIE,EBCMJK->aBCIJK');

D11 = -einsum_kg(H2B.ovvo(hA,PB,pA,HB),T3B.ppPhHH,'mBeJ,aeCmIK->aBCIJK')...
+einsum_kg(H2B.ovvo(hA,PB,PA,HB),T3B.PpPhHH,'mBEJ,EaCmIK->aBCIJK')...
+einsum_kg(H2B.ovvo(HA,PB,pA,HB),T3B.ppPHHH,'MBeJ,aeCIMK->aBCIJK')...
-einsum_kg(H2B.ovvo(HA,PB,PA,HB),T3B.PpPHHH,'MBEJ,EaCIMK->aBCIJK');
D11 = D11 - permute(D11,[1,3,2,4,5,6]) - permute(D11,[1,2,3,4,6,5]) + permute(D11,[1,3,2,4,6,5]);

D12 = -einsum_kg(H2C.voov(PB,hB,HB,pB),T3C.pPpHhH,'BmJe,aCeImK->aBCIJK')...
+einsum_kg(H2C.voov(PB,hB,HB,PB),T3C.pPPHhH,'BmJE,aECImK->aBCIJK')...
-einsum_kg(H2C.voov(PB,HB,HB,pB),T3C.pPpHHH,'BMJe,aCeIMK->aBCIJK')...
+einsum_kg(H2C.voov(PB,HB,HB,PB),T3C.pPPHHH,'BMJE,aECIMK->aBCIJK');
D12 = D12 - permute(D12,[1,3,2,4,5,6]) - permute(D12,[1,2,3,4,6,5]) + permute(D12,[1,3,2,4,6,5]);

D13 = +einsum_kg(H2B.ovov(hA,PB,HA,pB),T3C.pPphHH,'mBIe,aCemJK->aBCIJK')...
-einsum_kg(H2B.ovov(hA,PB,HA,PB),T3C.pPPhHH,'mBIE,aECmJK->aBCIJK')...
+einsum_kg(H2B.ovov(HA,PB,HA,pB),T3C.pPpHHH,'MBIe,aCeMJK->aBCIJK')...
-einsum_kg(H2B.ovov(HA,PB,HA,PB),T3C.pPPHHH,'MBIE,aECMJK->aBCIJK');
D13 = D13 - permute(D13,[1,3,2,4,5,6]);

D14 = -einsum_kg(H2B.vovo(pA,hB,pA,HB),T3C.pPPHhH,'ameJ,eBCImK->aBCIJK')...
-einsum_kg(H2B.vovo(pA,hB,PA,HB),T3C.PPPHhH,'amEJ,EBCImK->aBCIJK')...
-einsum_kg(H2B.vovo(pA,HB,pA,HB),T3C.pPPHHH,'aMeJ,eBCIMK->aBCIJK')...
-einsum_kg(H2B.vovo(pA,HB,PA,HB),T3C.PPPHHH,'aMEJ,EBCIMK->aBCIJK');
D14 = D14 - permute(D14,[1,2,3,4,6,5]);

X3C_pPPHHH = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14 + MM23C;



end

function [H1A,H1B,H2A,H2B,H2C,Vt3B,Vt3C] = get_t3c_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys)

    % h1A(mi)
    H1A.oo = sys.fa_oo ...
             +einsum_kg(sys.fa_ov,t1a,'me,ei->mi')...
             +einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi')...
             +einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi')...
             +0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi')...
             +einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi')...
             +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t1a,'me,ei->mi')...
             +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t1a,'me,ei->mi');
         
    % h1A(ae)
    H1A.vv = sys.fa_vv ...
             -einsum_kg(sys.fa_ov,t1a,'me,am->ae')...
             +einsum_kg(sys.vA_vovv,t1a,'amef,fm->ae')...
             +einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae')...
             -0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,afmn->ae')...
             -einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ae')...
             -einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t1a,'me,am->ae')...
             -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t1a,'me,am->ae');
         
    % h1A(me)
    H1A.ov = sys.fa_ov + einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me') + einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');

    % h1B(mi)
    H1B.oo = sys.fb_oo ...
             +einsum_kg(sys.fb_ov,t1b,'me,ei->mi')...
             +einsum_kg(sys.vB_oovo,t1a,'nmfi,fn->mi')...
             +einsum_kg(sys.vC_ooov,t1b,'mnif,fn->mi')...
             +einsum_kg(sys.vB_oovv,t2b,'nmfe,feni->mi')...
             +0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,efin->mi')...
             +einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t1b,'me,ei->mi')...
             +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t1b,'me,ei->mi');
         
    % h1B(ae)
    H1B.vv = sys.fb_vv ...
             -einsum_kg(sys.fb_ov,t1b,'me,am->ae')...
             +einsum_kg(sys.vB_ovvv,t1a,'nafe,fn->ae')...
             +einsum_kg(sys.vC_vovv,t1b,'anef,fn->ae')...
             -einsum_kg(sys.vB_oovv,t2b,'nmfe,fanm->ae')...
             -0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,afmn->ae')...
             -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t1b,'me,am->ae')...
             -einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t1b,'me,am->ae');
         
    % h1B(me)
    H1B.ov = sys.fb_ov + einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me') + einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
    
    % h2A(amie)
    H2A.voov = sys.vA_voov ...
                -einsum_kg(sys.vA_ooov,t1a,'nmie,an->amie')...
                +einsum_kg(sys.vA_vovv,t1a,'amfe,fi->amie')...
                -einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie')...
                +einsum_kg(sys.vA_oovv,t2a,'nmfe,afin->amie')...
                +einsum_kg(sys.vB_oovv,t2b,'mnef,afin->amie');

     % h2B(amij)
     H2B.vooo = sys.vB_vooo...
               -einsum_kg(sys.vB_oooo,t1a,'nmij,an->amij')...
               +einsum_kg(sys.vB_vovo,t1a,'amej,ei->amij')...
               +einsum_kg(sys.vB_voov,t1b,'amie,ej->amij')...
               +einsum_kg(sys.vB_oovo,t2a,'nmej,aein->amij')...
               +einsum_kg(sys.vC_oovo,t2b,'nmej,aein->amij')...
               +einsum_kg(sys.vB_vovv,t2b,'amef,efij->amij')...
               -einsum_kg(sys.vB_ooov,t2b,'nmif,afnj->amij')...
               +einsum_kg(sys.fb_ov,t2b,'me,aeij->amij')...
               -einsum_kg(einsum_kg(sys.vB_oovo,t1a,'nmej,an->amej'),t1a,'amej,ei->amij')...
               -einsum_kg(einsum_kg(sys.vB_ooov,t1b,'nmie,ej->nmij'),t1a,'nmij,an->amij')...
               +einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ej->mnjf'),t2b,'mnjf,afin->amij')...
               +einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t2b,'me,aeij->amij')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,an->amfe'),t2b,'amfe,feij->amij')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t2b,'nmie,aenj->amij')...
               +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t2b,'me,aeij->amij')...
               +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,ej->nmfj'),t2a,'nmfj,afin->amij')...
               -einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,ej->nmfj'),t1a,'nmfj,an->amfj'),t1a,'amfj,fi->amij')...
               +einsum_kg(einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie'),t1b,'amie,ej->amij');

    % h2B(maji)
    H2B.ovoo =  sys.vB_ovoo...
                +einsum_kg(sys.vB_ovvo,t1a,'maei,ej->maji')...
                -einsum_kg(sys.vB_oooo,t1b,'mnji,an->maji')...
                +einsum_kg(sys.vB_ovov,t1b,'majf,fi->maji')...
                +einsum_kg(sys.vB_ooov,t2c,'mnjf,afin->maji')...
                +einsum_kg(sys.vB_ovvv,t2b,'maef,efji->maji')...
                +einsum_kg(sys.vA_ooov,t2b,'mnjf,fani->maji')...
                -einsum_kg(sys.vB_oovo,t2b,'mnei,eajn->maji')...
                +einsum_kg(sys.fa_ov,t2b,'me,eaji->maji')...
                -einsum_kg(einsum_kg(sys.vB_ooov,t1b,'mnjf,an->majf'),t1b,'majf,fi->maji')...
                +einsum_kg(einsum_kg(sys.vB_ovvv,t1a,'maef,ej->majf'),t1b,'majf,fi->maji')...
                -einsum_kg(einsum_kg(sys.vB_oovo,t1a,'mnei,ej->mnji'),t1b,'mnji,an->maji')...
                +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ej->mnjf'),t2c,'mnjf,afin->maji')...
                +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t2b,'me,eaji->maji')...
                +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ej->mnjf'),t2b,'mnjf,fani->maji')...
                -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,an->maef'),t2b,'maef,efji->maji')...
                +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t2b,'me,eaji->maji')...
                -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t2b,'mnei,eajn->maji')...
                -einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ej->mnjf'),t1b,'mnjf,an->majf'),t1b,'majf,fi->maji');

    % h2B(abie)
    H2B.vvov =  sys.vB_vvov...
	    	-einsum_kg(sys.vB_ovov,t1a,'nbie,an->abie')...
	    	+einsum_kg(sys.vB_vvvv,t1a,'abfe,fi->abie')...
	    	-einsum_kg(sys.vB_voov,t1b,'amie,bm->abie')...
	    	+einsum_kg(sys.vB_ovvv,t2a,'nbfe,afin->abie')...
	    	+einsum_kg(sys.vC_ovvv,t2b,'nbfe,afin->abie')...
	    	+einsum_kg(sys.vB_ooov,t2b,'nmie,abnm->abie')...
	    	-einsum_kg(sys.vB_vovv,t2b,'amfe,fbim->abie')...
	    	-einsum_kg(einsum_kg(sys.vB_ovvv,t1a,'nbfe,an->abfe'),t1a,'abfe,fi->abie')...
	    	+einsum_kg(einsum_kg(sys.vB_ooov,t1a,'nmie,an->amie'),t1b,'amie,bm->abie')...
	    	-einsum_kg(einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie'),t1b,'amie,bm->abie')...
	    	-einsum_kg(sys.fb_ov,t2b,'me,abim->abie')...
	    	-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bm->bnef'),t2b,'bnef,afin->abie')...
	    	-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t2b,'me,abim->abie')...
	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t2b,'nmie,abnm->abie')...
	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,an->amfe'),t2b,'amfe,fbim->abie')...
	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t2b,'me,abim->abie')...
	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,bm->nbfe'),t2a,'nbfe,afin->abie')...
	    	+einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,an->amfe'),t1a,'amfe,fi->amie'),t1b,'amie,bm->abie');

    % h2B(baei)
    H2B.vvvo =  sys.vB_vvvo...
	    	+einsum_kg(sys.vB_vvvv,t1b,'baef,fi->baei')...
	    	-einsum_kg(sys.vB_vovo,t1b,'bnei,an->baei')...
	    	-einsum_kg(sys.vB_ovvo,t1a,'maei,bm->baei')...
	    	+einsum_kg(sys.vB_vovv,t2c,'bnef,afin->baei')...
	    	+einsum_kg(sys.vB_oovo,t2b,'mnei,bamn->baei')...
	    	+einsum_kg(sys.vA_vovv,t2b,'bnef,fani->baei')...
	    	-einsum_kg(sys.vB_ovvv,t2b,'maef,bfmi->baei')...
	    	-einsum_kg(sys.fa_ov,t2b,'me,bami->baei')...
	    	-einsum_kg(einsum_kg(sys.vB_vovv,t1b,'bnef,fi->bnei'),t1b,'bnei,an->baei')...
	    	+einsum_kg(einsum_kg(sys.vB_oovo,t1b,'mnei,an->maei'),t1a,'maei,bm->baei')...
	    	-einsum_kg(einsum_kg(sys.vB_ovvv,t1b,'maef,fi->maei'),t1a,'maei,bm->baei')...
	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,bm->bnef'),t2c,'bnef,afin->baei')...
	    	-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t2b,'me,bami->baei')...
	    	-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bm->bnef'),t2b,'bnef,fani->baei')...
	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t2b,'mnei,bamn->baei')...
	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t2b,'me,bami->baei')...
	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,an->maef'),t2b,'maef,bfmi->baei')...
	    	+einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,an->maef'),t1b,'maef,fi->maei'),t1a,'maei,bm->baei');

    % h2B(mnij)
    H2B.oooo =  sys.vB_oooo...
	    	+einsum_kg(sys.vB_oovo,t1a,'mnej,ei->mnij')...
	   	+einsum_kg(sys.vB_ooov,t1b,'mnif,fj->mnij')...
	   	+einsum_kg(sys.vB_oovv,t2b,'mnef,efij->mnij')...
	   	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ei->mnif'),t1b,'mnif,fj->mnij'); 

    % h2B(abef)
    H2B.vvvv =  sys.vB_vvvv...
	    	-einsum_kg(sys.vB_ovvv,t1a,'mbef,am->abef')...
	    	-einsum_kg(sys.vB_vovv,t1b,'anef,bn->abef')...
	    	+einsum_kg(sys.vB_oovv,t2b,'mnef,abmn->abef')...
	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,am->anef'),t1b,'anef,bn->abef');
        
    % h2B(amie)
    H2B.voov = sys.vB_voov ...
               -einsum_kg(sys.vB_ooov,t1a,'nmie,an->amie')...
               +einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie')...
               +einsum_kg(sys.vB_oovv,t2a,'nmfe,afin->amie')...
               +einsum_kg(sys.vC_oovv,t2b,'mnef,afin->amie');
          
    % h2B(maie)
    H2B.ovov = sys.vB_ovov ...
               +einsum_kg(sys.vB_ovvv,t1a,'mafe,fi->maie')...
               -einsum_kg(sys.vB_ooov,t1b,'mnie,an->maie')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnfe,an->mafe'),t1a,'mafe,fi->maie')...
               -einsum_kg(sys.vB_oovv,t2b,'mnfe,fain->maie');
           
    % h2B(amei)
    H2B.vovo = sys.vB_vovo ...
               -einsum_kg(sys.vB_oovo,t1a,'nmei,an->amei')...
               +einsum_kg(sys.vB_vovv,t1b,'amef,fi->amei')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmef,fi->nmei'),t1a,'nmei,an->amei')...
               -einsum_kg(sys.vB_oovv,t2b,'nmef,afni->amei');
           
    % h2B(maei)
    H2B.ovvo = sys.vB_ovvo ...
               +einsum_kg(sys.vB_ovvv,t1b,'maef,fi->maei')...
               -einsum_kg(sys.vB_oovo,t1b,'mnei,an->maei')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->maei')...
               +einsum_kg(sys.vB_oovv,t2c,'mnef,afin->maei')...
               +einsum_kg(sys.vA_oovv,t2b,'mnef,fani->maei');
           
    % h2C(amij)
    H2C.vooo =  sys.vC_vooo ...
	    	+einsum_kg(sys.vC_vovo,t1b,'amej,ei->amij')...
	    			-einsum_kg(sys.vC_vovo,t1b,'amei,ej->amij')...
	    	-einsum_kg(sys.vC_oooo,t1b,'nmij,an->amij')...
	    	+einsum_kg(sys.vB_oovo,t2b,'nmfj,fani->amij')...
	    			-einsum_kg(sys.vB_oovo,t2b,'nmfi,fanj->amij')...
	    	+0.5*einsum_kg(sys.vC_vovv,t2c,'amfe,feij->amij')...
	    	+einsum_kg(sys.vC_ooov,t2c,'mnjf,afin->amij')...
	    			-einsum_kg(sys.vC_ooov,t2c,'mnif,afjn->amij')...
	    	+einsum_kg(sys.fb_ov,t2c,'me,aeij->amij')...
	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t2c,'me,aeij->amij')...
	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,ej->nmfj'),t2b,'nmfj,fani->amij')...
	    			-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,ei->nmfi'),t2b,'nmfi,fanj->amij')...
	    	-einsum_kg(einsum_kg(sys.vC_ooov,t1b,'mnjf,an->majf'),t1b,'majf,fi->amij')...
	    			+einsum_kg(einsum_kg(sys.vC_ooov,t1b,'mnif,an->maif'),t1b,'maif,fj->amij')...
	    	+einsum_kg(einsum_kg(sys.vC_vovv,t1b,'amfe,fi->amie'),t1b,'amie,ej->amij')...
	    	-0.5*einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,an->maef'),t2c,'maef,feij->amij')...
	    	+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ej->mnjf'),t2c,'mnjf,afin->amij')...
	    			-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ei->mnif'),t2c,'mnif,afjn->amij')...
	    	+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t2c,'me,aeij->amij')...
	    	-einsum_kg(einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ej->mnjf'),t1b,'mnjf,fi->mnji'),t1b,'mnji,an->amij');
    H2C.ovoo = -permute(H2C.vooo,[2,1,3,4]);
           
    % h2C(abie)
    H2C.vvov =  sys.vC_vvov...
	    	+einsum_kg(sys.vC_vvvv,t1b,'abfe,fi->abie')...
	    	-einsum_kg(sys.vC_ovov,t1b,'nbie,an->abie')...
	    			+einsum_kg(sys.vC_ovov,t1b,'naie,bn->abie')...
	    	+einsum_kg(sys.vB_ovvv,t2b,'nbfe,fani->abie')...
	    			-einsum_kg(sys.vB_ovvv,t2b,'nafe,fbni->abie')...
	    	+0.5*einsum_kg(sys.vC_oovo,t2c,'mnei,abnm->abie')...
	    	+einsum_kg(sys.vC_ovvv,t2c,'nbfe,afin->abie')...
	    			-einsum_kg(sys.vC_ovvv,t2c,'nafe,bfin->abie')...
	    	-einsum_kg(sys.fb_ov,t2c,'me,abim->abie')...
	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t2c,'me,abim->abie')...
	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,bm->nbfe'),t2b,'nbfe,fani->abie')...
	    			+einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,am->nafe'),t2b,'nafe,fbni->abie')...
	    	-einsum_kg(einsum_kg(sys.vC_ovvv,t1b,'nbfe,an->abfe'),t1b,'abfe,fi->abie')...
	    			+einsum_kg(einsum_kg(sys.vC_ovvv,t1b,'nafe,bn->bafe'),t1b,'bafe,fi->abie')...
	    	+einsum_kg(einsum_kg(sys.vC_oovo,t1b,'mnei,bm->bnei'),t1b,'bnei,an->abie')...
	    	-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bm->bnef'),t2c,'bnef,afin->abie')...
	    			+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,am->anef'),t2c,'anef,bfin->abie')...
	    	+0.5*einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fi->mnei'),t2c,'mnei,abnm->abie')...
	    	-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t2c,'me,abim->abie')...
	    	+einsum_kg(einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bm->bnef'),t1b,'bnef,an->baef'),t1b,'baef,fi->abie');
     H2C.vvvo = -permute(H2C.vvov,[1,2,4,3]);       
     
     % h2C(mnij)
     H2C.oooo = sys.vC_oooo...
	      	+einsum_kg(sys.vC_ooov,t1b,'mnie,ej->mnij')...
	    			-einsum_kg(sys.vC_ooov,t1b,'mnje,ei->mnij')...
	    	+0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,efij->mnij')...
	    	+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ei->mnif'),t1b,'mnif,fj->mnij'); 

     % h2C(abef)
     H2C.vvvv = sys.vC_vvvv...
	     	-einsum_kg(sys.vC_ovvv,t1b,'mbef,am->abef')...
	     			+einsum_kg(sys.vC_ovvv,t1b,'maef,bm->abef')...
	     	+0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,abmn->abef')...
	     	+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bn->mbef'),t1b,'mbef,am->abef');
	      
     % h2C(amie)
     H2C.voov = sys.vC_voov ...
                -einsum_kg(sys.vC_oovo,t1b,'mnei,an->amie')...
                +einsum_kg(sys.vC_vovv,t1b,'amfe,fi->amie')...
                -einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->amie')...
                +einsum_kg(sys.vC_oovv,t2c,'mnef,afin->amie')...
                +einsum_kg(sys.vB_oovv,t2b,'nmfe,fani->amie');

     % (Vt3)_C intermediates
     PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
     PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
     
d1 = -einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3B.ppPhHh,'nmfe,afBnIm->aBIe')...
+einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3B.PpPhHh,'nmFe,FaBnIm->aBIe')...
+einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3B.ppPHHh,'Nmfe,afBINm->aBIe')...
-einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3B.PpPHHh,'NmFe,FaBINm->aBIe')...
-einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3B.ppPhHH,'nMfe,afBnIM->aBIe')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3B.PpPhHH,'nMFe,FaBnIM->aBIe')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3B.ppPHHH,'NMfe,afBINM->aBIe')...
-einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3B.PpPHHH,'NMFe,FaBINM->aBIe');

d2 = -0.5*einsum_kg(sys.vC_oovv(hB,hB,pB,:),T3C.pPpHhh,'nmfe,aBfInm->aBIe')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,PB,:),T3C.pPPHhh,'nmFe,aFBInm->aBIe')...
-einsum_kg(sys.vC_oovv(hB,HB,pB,:),T3C.pPpHhH,'nMfe,aBfInM->aBIe')...
+einsum_kg(sys.vC_oovv(hB,HB,PB,:),T3C.pPPHhH,'nMFe,aFBInM->aBIe')...
-0.5*einsum_kg(sys.vC_oovv(HB,HB,pB,:),T3C.pPpHHH,'NMfe,aBfINM->aBIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),T3C.pPPHHH,'NMFe,aFBINM->aBIe');

Vt3B.pPHv = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3B.ppphHH,'nmfe,afenIJ->amIJ')...
-einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3B.ppPhHH,'nmfE,afEnIJ->amIJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3B.PpphHH,'nmFe,FaenIJ->amIJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3B.PpPhHH,'nmFE,FaEnIJ->amIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3B.pppHHH,'Nmfe,afeINJ->amIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3B.ppPHHH,'NmfE,afEINJ->amIJ')...
-einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3B.PppHHH,'NmFe,FaeINJ->amIJ')...
-einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3B.PpPHHH,'NmFE,FaEINJ->amIJ');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,:,pB,pB),T3C.pppHhH,'nmfe,afeInJ->amIJ')...
+einsum_kg(sys.vC_oovv(hB,:,PB,pB),T3C.pPpHhH,'nmFe,aFeInJ->amIJ')...
+0.5*einsum_kg(sys.vC_oovv(hB,:,PB,PB),T3C.pPPHhH,'nmFE,aFEInJ->amIJ')...
+0.5*einsum_kg(sys.vC_oovv(HB,:,pB,pB),T3C.pppHHH,'Nmfe,afeINJ->amIJ')...
+einsum_kg(sys.vC_oovv(HB,:,PB,pB),T3C.pPpHHH,'NmFe,aFeINJ->amIJ')...
+0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),T3C.pPPHHH,'NmFE,aFEINJ->amIJ');

Vt3B.poHH = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3C.pPPhhH,'nmfe,fABnmI->BAeI')...
-einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3C.PPPhhH,'nmFe,FABnmI->BAeI')...
-einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3C.pPPHhH,'Nmfe,fABNmI->BAeI')...
-einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3C.PPPHhH,'NmFe,FABNmI->BAeI')...
+einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3C.pPPhHH,'nMfe,fABnIM->BAeI')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3C.PPPhHH,'nMFe,FABnIM->BAeI')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3C.pPPHHH,'NMfe,fABNIM->BAeI')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPPHHH,'NMFe,FABNIM->BAeI');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,hB,:,pB),T3D.PPphhH,'mnef,ABfmnI->BAeI')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,:,PB),T3D.PPPhhH,'mneF,ABFmnI->BAeI')...
+einsum_kg(sys.vC_oovv(HB,hB,:,pB),T3D.PPphHH,'Mnef,ABfnIM->BAeI')...
+einsum_kg(sys.vC_oovv(HB,hB,:,PB),T3D.PPPhHH,'MneF,ABFnIM->BAeI')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,pB),T3D.PPpHHH,'MNef,ABfIMN->BAeI')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPPHHH,'MNeF,ABFIMN->BAeI');

Vt3C.PPvH = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3C.pPphHH,'nmfe,fAenIJ->mAJI')...
+einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3C.pPPhHH,'nmfE,fAEnIJ->mAJI')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3C.PPphHH,'nmFe,FAenIJ->mAJI')...
+einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3C.PPPhHH,'nmFE,FAEnIJ->mAJI')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3C.pPpHHH,'Nmfe,fAeNIJ->mAJI')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3C.pPPHHH,'NmfE,fAENIJ->mAJI')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3C.PPpHHH,'NmFe,FAeNIJ->mAJI')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPPHHH,'NmFE,FAENIJ->mAJI');

d2 = +0.5*einsum_kg(sys.vC_oovv(:,hB,pB,pB),T3D.PpphHH,'mnef,AefnIJ->mAJI')...
-einsum_kg(sys.vC_oovv(:,hB,pB,PB),T3D.PPphHH,'mneF,AFenIJ->mAJI')...
+0.5*einsum_kg(sys.vC_oovv(:,hB,PB,PB),T3D.PPPhHH,'mnEF,AEFnIJ->mAJI')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,pB,pB),T3D.PppHHH,'mNef,AefIJN->mAJI')...
-einsum_kg(sys.vC_oovv(:,HB,pB,PB),T3D.PPpHHH,'mNeF,AFeIJN->mAJI')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPPHHH,'mNEF,AEFIJN->mAJI');

Vt3C.oPHH = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3C.pPphhH,'mnef,aBfmnJ->aBeJ')...
+einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3C.pPPhhH,'mneF,aFBmnJ->aBeJ')...
-einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3C.pPphHH,'mNef,aBfmNJ->aBeJ')...
+einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3C.pPPhHH,'mNeF,aFBmNJ->aBeJ')...
-einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3C.pPpHhH,'Mnef,aBfMnJ->aBeJ')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3C.pPPHhH,'MneF,aFBMnJ->aBeJ')...
-einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3C.pPpHHH,'MNef,aBfMNJ->aBeJ')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.pPPHHH,'MNeF,aFBMNJ->aBeJ');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3B.ppPhhH,'mnef,afBmnJ->aBeJ')...
-0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3B.PpPhhH,'mneF,FaBmnJ->aBeJ')...
-einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3B.ppPhHH,'Mnef,afBnMJ->aBeJ')...
+einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3B.PpPhHH,'MneF,FaBnMJ->aBeJ')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3B.ppPHHH,'MNef,afBMNJ->aBeJ')...
-0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PpPHHH,'MNeF,FaBMNJ->aBeJ');

Vt3B.pPvH = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3C.pPpHhH,'mnef,eBfInJ->mBIJ')...
-einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3C.PPpHhH,'mnEf,EBfInJ->mBIJ')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3C.pPPHhH,'mneF,eFBInJ->mBIJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3C.PPPHhH,'mnEF,EFBInJ->mBIJ')...
-einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3C.pPpHHH,'mNef,eBfINJ->mBIJ')...
-einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PPpHHH,'mNEf,EBfINJ->mBIJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPPHHH,'mNeF,eFBINJ->mBIJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPPHHH,'mNEF,EFBINJ->mBIJ');

d2 = -0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3B.ppPhHH,'mnef,efBnIJ->mBIJ')...
+einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3B.PpPhHH,'mneF,FeBnIJ->mBIJ')...
-0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3B.PPPhHH,'mnEF,EFBnIJ->mBIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3B.ppPHHH,'mNef,efBINJ->mBIJ')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PpPHHH,'mNeF,FeBINJ->mBIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPPHHH,'mNEF,EFBINJ->mBIJ');

Vt3B.oPHH = d1 + d2;

     
end