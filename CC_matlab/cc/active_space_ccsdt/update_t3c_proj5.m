function [X3C_PPPhHH] = update_t3c_proj5(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, sys, shift)
    
    [H1A,H1B,H2A,H2B,H2C,Vt3B,Vt3C] = get_t3c_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys);
    
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
    
    % MM23C
    I2B_vooo = H2B.vooo - einsum_kg(H1B.ov,t2b,'me,aeij->amij');
    I2C_ovoo = H2C.ovoo - einsum_kg(H1B.ov,t2c,'me,ecjk->mcjk');
    I2B_ovoo = H2B.ovoo - einsum_kg(H1A.ov,t2b,'me,ebij->mbij');
    

    M23_D1 = +einsum_kg(H2B.vvov(PA,PB,hA,:) + Vt3B.PPhv,t2c(:,PB,HB,HB),'ABie,eCJK->ABCiJK');
    
    M23_D2 = -einsum_kg(I2B_vooo(PA,:,hA,HB) + Vt3B.PohH,t2c(PB,PB,:,HB),'AmiJ,BCmK->ABCiJK');
  
    M23_D3 = +einsum_kg(H2C.vvvo(PB,PB,:,HB) + Vt3C.PPvH,t2b(PA,:,hA,HB),'BCeK,AeiJ->ABCiJK');

    M23_D4 = -einsum_kg(I2C_ovoo(:,PB,HB,HB) + Vt3C.oPHH,t2b(PA,PB,hA,:),'mCJK,ABim->ABCiJK');

    M23_D5 = +einsum_kg(H2B.vvvo(PA,PB,:,HB) + Vt3B.PPvH,t2b(:,PB,hA,HB),'ABeJ,eCiK->ABCiJK');

    M23_D6 = -einsum_kg(I2B_ovoo(:,PB,hA,HB) + Vt3B.oPhH,t2b(PA,PB,:,HB),'mBiJ,ACmK->ABCiJK');

    M23_14 = M23_D1 + M23_D4;
    M23_14 = M23_14 - permute(M23_14,[1,3,2,4,5,6]);
    M23_23 = M23_D2 + M23_D3;
    M23_23 = M23_23 - permute(M23_23,[1,2,3,4,6,5]);
    M23_56 = M23_D5 + M23_D6;
    M23_56 = M23_56 - permute(M23_56,[1,3,2,4,5,6]) - permute(M23_56,[1,2,3,4,6,5]) + permute(M23_56,[1,3,2,4,6,5]);
    MM23C = M23_14 + M23_23 + M23_56;

    % (HBar*T3)_C
D1 = -einsum_kg(H1A.oo(hA,hA),T3C.PPPhHH,'mi,ABCmJK->ABCiJK')...
-einsum_kg(H1A.oo(HA,hA),T3C.PPPHHH,'Mi,ABCMJK->ABCiJK');

D2 = -einsum_kg(H1B.oo(hB,HB),T3C.PPPhhH,'mJ,ABCimK->ABCiJK')...
-einsum_kg(H1B.oo(HB,HB),T3C.PPPhHH,'MJ,ABCiMK->ABCiJK');
D2 = D2 - permute(D2,[1,2,3,4,6,5]);

D3 = +einsum_kg(H1A.vv(PA,pA),T3C.pPPhHH,'Ae,eBCiJK->ABCiJK')...
+einsum_kg(H1A.vv(PA,PA),T3C.PPPhHH,'AE,EBCiJK->ABCiJK');

D4 = -einsum_kg(H1B.vv(PB,pB),T3C.PPphHH,'Be,ACeiJK->ABCiJK')...
+einsum_kg(H1B.vv(PB,PB),T3C.PPPhHH,'BE,AECiJK->ABCiJK');
D4 = D4 - permute(D4,[1,3,2,4,5,6]);

D5 = +0.5*einsum_kg(H2C.oooo(hB,hB,HB,HB),T3C.PPPhhh,'mnJK,ABCimn->ABCiJK')...
-einsum_kg(H2C.oooo(HB,hB,HB,HB),T3C.PPPhhH,'MnJK,ABCinM->ABCiJK')...
+0.5*einsum_kg(H2C.oooo(HB,HB,HB,HB),T3C.PPPhHH,'MNJK,ABCiMN->ABCiJK');

D6 = +einsum_kg(H2B.oooo(hA,hB,hA,HB),T3C.PPPhhH,'mniJ,ABCmnK->ABCiJK')...
+einsum_kg(H2B.oooo(hA,HB,hA,HB),T3C.PPPhHH,'mNiJ,ABCmNK->ABCiJK')...
+einsum_kg(H2B.oooo(HA,hB,hA,HB),T3C.PPPHhH,'MniJ,ABCMnK->ABCiJK')...
+einsum_kg(H2B.oooo(HA,HB,hA,HB),T3C.PPPHHH,'MNiJ,ABCMNK->ABCiJK');
D6 = D6 - permute(D6,[1,2,3,4,6,5]);

D7 = +0.5*einsum_kg(H2C.vvvv(PB,PB,pB,pB),T3C.PpphHH,'BCef,AefiJK->ABCiJK')...
+einsum_kg(H2C.vvvv(PB,PB,PB,pB),T3C.PPphHH,'BCEf,AEfiJK->ABCiJK')...
+0.5*einsum_kg(H2C.vvvv(PB,PB,PB,PB),T3C.PPPhHH,'BCEF,AEFiJK->ABCiJK');

D8 = -einsum_kg(H2B.vvvv(PA,PB,pA,pB),T3C.pPphHH,'ABef,eCfiJK->ABCiJK')...
+einsum_kg(H2B.vvvv(PA,PB,pA,PB),T3C.pPPhHH,'ABeF,eFCiJK->ABCiJK')...
-einsum_kg(H2B.vvvv(PA,PB,PA,pB),T3C.PPphHH,'ABEf,ECfiJK->ABCiJK')...
+einsum_kg(H2B.vvvv(PA,PB,PA,PB),T3C.PPPhHH,'ABEF,EFCiJK->ABCiJK');
D8 = D8 - permute(D8,[1,3,2,4,5,6]);

D9 = +einsum_kg(H2A.voov(PA,hA,hA,pA),T3C.pPPhHH,'Amie,eBCmJK->ABCiJK')...
+einsum_kg(H2A.voov(PA,hA,hA,PA),T3C.PPPhHH,'AmiE,EBCmJK->ABCiJK')...
+einsum_kg(H2A.voov(PA,HA,hA,pA),T3C.pPPHHH,'AMie,eBCMJK->ABCiJK')...
+einsum_kg(H2A.voov(PA,HA,hA,PA),T3C.PPPHHH,'AMiE,EBCMJK->ABCiJK');

D10 = +einsum_kg(H2B.voov(PA,hB,hA,pB),T3D.PPphHH,'Amie,BCemJK->ABCiJK')...
+einsum_kg(H2B.voov(PA,hB,hA,PB),T3D.PPPhHH,'AmiE,EBCmJK->ABCiJK')...
+einsum_kg(H2B.voov(PA,HB,hA,pB),T3D.PPpHHH,'AMie,BCeMJK->ABCiJK')...
+einsum_kg(H2B.voov(PA,HB,hA,PB),T3D.PPPHHH,'AMiE,EBCMJK->ABCiJK');

D11 = +einsum_kg(H2B.ovvo(hA,PB,pA,HB),T3B.PpPhhH,'mBeJ,AeCimK->ABCiJK')...
+einsum_kg(H2B.ovvo(hA,PB,PA,HB),T3B.PPPhhH,'mBEJ,AECimK->ABCiJK')...
+einsum_kg(H2B.ovvo(HA,PB,pA,HB),T3B.PpPhHH,'MBeJ,AeCiMK->ABCiJK')...
+einsum_kg(H2B.ovvo(HA,PB,PA,HB),T3B.PPPhHH,'MBEJ,AECiMK->ABCiJK');
D11 = D11 - permute(D11,[1,3,2,4,5,6]) - permute(D11,[1,2,3,4,6,5]) + permute(D11,[1,3,2,4,6,5]);

D12 = -einsum_kg(H2C.voov(PB,hB,HB,pB),T3C.PPphhH,'BmJe,ACeimK->ABCiJK')...
+einsum_kg(H2C.voov(PB,hB,HB,PB),T3C.PPPhhH,'BmJE,AECimK->ABCiJK')...
-einsum_kg(H2C.voov(PB,HB,HB,pB),T3C.PPphHH,'BMJe,ACeiMK->ABCiJK')...
+einsum_kg(H2C.voov(PB,HB,HB,PB),T3C.PPPhHH,'BMJE,AECiMK->ABCiJK');
D12 = D12 - permute(D12,[1,3,2,4,5,6]) - permute(D12,[1,2,3,4,6,5]) + permute(D12,[1,3,2,4,6,5]);

D13 = +einsum_kg(H2B.ovov(hA,PB,hA,pB),T3C.PPphHH,'mBie,ACemJK->ABCiJK')...
-einsum_kg(H2B.ovov(hA,PB,hA,PB),T3C.PPPhHH,'mBiE,AECmJK->ABCiJK')...
+einsum_kg(H2B.ovov(HA,PB,hA,pB),T3C.PPpHHH,'MBie,ACeMJK->ABCiJK')...
-einsum_kg(H2B.ovov(HA,PB,hA,PB),T3C.PPPHHH,'MBiE,AECMJK->ABCiJK');
D13 = D13 - permute(D13,[1,3,2,4,5,6]);

D14 = -einsum_kg(H2B.vovo(PA,hB,pA,HB),T3C.pPPhhH,'AmeJ,eBCimK->ABCiJK')...
-einsum_kg(H2B.vovo(PA,hB,PA,HB),T3C.PPPhhH,'AmEJ,EBCimK->ABCiJK')...
-einsum_kg(H2B.vovo(PA,HB,pA,HB),T3C.pPPhHH,'AMeJ,eBCiMK->ABCiJK')...
-einsum_kg(H2B.vovo(PA,HB,PA,HB),T3C.PPPhHH,'AMEJ,EBCiMK->ABCiJK');
D14 = D14 - permute(D14,[1,2,3,4,6,5]);

X3C_PPPhHH = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14 + MM23C;



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
     
d1 = +einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3B.PpPhhh,'nmfe,AfBinm->ABie')...
+einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3B.PPPhhh,'nmFe,AFBinm->ABie')...
+einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3B.PpPhHh,'Nmfe,AfBiNm->ABie')...
+einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3B.PPPhHh,'NmFe,AFBiNm->ABie')...
+einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3B.PpPhhH,'nMfe,AfBinM->ABie')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3B.PPPhhH,'nMFe,AFBinM->ABie')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3B.PpPhHH,'NMfe,AfBiNM->ABie')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3B.PPPhHH,'NMFe,AFBiNM->ABie');

d2 = -0.5*einsum_kg(sys.vC_oovv(hB,hB,pB,:),T3C.PPphhh,'nmfe,ABfinm->ABie')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,PB,:),T3C.PPPhhh,'nmFe,AFBinm->ABie')...
-einsum_kg(sys.vC_oovv(hB,HB,pB,:),T3C.PPphhH,'nMfe,ABfinM->ABie')...
+einsum_kg(sys.vC_oovv(hB,HB,PB,:),T3C.PPPhhH,'nMFe,AFBinM->ABie')...
-0.5*einsum_kg(sys.vC_oovv(HB,HB,pB,:),T3C.PPphHH,'NMfe,ABfiNM->ABie')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),T3C.PPPhHH,'NMFe,AFBiNM->ABie');

Vt3B.PPhv = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3B.PpphhH,'nmfe,AfeinJ->AmiJ')...
+einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3B.PpPhhH,'nmfE,AfEinJ->AmiJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3B.PPphhH,'nmFe,AFeinJ->AmiJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3B.PPPhhH,'nmFE,AFEinJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3B.PpphHH,'Nmfe,AfeiNJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3B.PpPhHH,'NmfE,AfEiNJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3B.PPphHH,'NmFe,AFeiNJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3B.PPPhHH,'NmFE,AFEiNJ->AmiJ');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,:,pB,pB),T3C.PpphhH,'nmfe,AfeinJ->AmiJ')...
+einsum_kg(sys.vC_oovv(hB,:,PB,pB),T3C.PPphhH,'nmFe,AFeinJ->AmiJ')...
+0.5*einsum_kg(sys.vC_oovv(hB,:,PB,PB),T3C.PPPhhH,'nmFE,AFEinJ->AmiJ')...
+0.5*einsum_kg(sys.vC_oovv(HB,:,pB,pB),T3C.PpphHH,'Nmfe,AfeiNJ->AmiJ')...
+einsum_kg(sys.vC_oovv(HB,:,PB,pB),T3C.PPphHH,'NmFe,AFeiNJ->AmiJ')...
+0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),T3C.PPPhHH,'NmFE,AFEiNJ->AmiJ');

Vt3B.PohH = d1 + d2;

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

d1 = -einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3C.PPphhH,'mnef,ABfmnJ->ABeJ')...
+einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3C.PPPhhH,'mneF,AFBmnJ->ABeJ')...
-einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3C.PPphHH,'mNef,ABfmNJ->ABeJ')...
+einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3C.PPPhHH,'mNeF,AFBmNJ->ABeJ')...
-einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3C.PPpHhH,'Mnef,ABfMnJ->ABeJ')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3C.PPPHhH,'MneF,AFBMnJ->ABeJ')...
-einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3C.PPpHHH,'MNef,ABfMNJ->ABeJ')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.PPPHHH,'MNeF,AFBMNJ->ABeJ');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3B.PpPhhH,'mnef,AfBmnJ->ABeJ')...
+0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3B.PPPhhH,'mneF,AFBmnJ->ABeJ')...
-einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3B.PpPhHH,'Mnef,AfBnMJ->ABeJ')...
-einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3B.PPPhHH,'MneF,AFBnMJ->ABeJ')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3B.PpPHHH,'MNef,AfBMNJ->ABeJ')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PPPHHH,'MNeF,AFBMNJ->ABeJ');

Vt3B.PPvH = -d1 - d2;
     

d1 = -einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3C.pPphhH,'mnef,eBfinJ->mBiJ')...
-einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3C.PPphhH,'mnEf,EBfinJ->mBiJ')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3C.pPPhhH,'mneF,eFBinJ->mBiJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3C.PPPhhH,'mnEF,EFBinJ->mBiJ')...
-einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3C.pPphHH,'mNef,eBfiNJ->mBiJ')...
-einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PPphHH,'mNEf,EBfiNJ->mBiJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPPhHH,'mNeF,eFBiNJ->mBiJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPPhHH,'mNEF,EFBiNJ->mBiJ');

d2 = +0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3B.ppPhhH,'mnef,efBinJ->mBiJ')...
-einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3B.PpPhhH,'mneF,FeBinJ->mBiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3B.PPPhhH,'mnEF,EFBinJ->mBiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3B.ppPhHH,'mNef,efBiNJ->mBiJ')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PpPhHH,'mNeF,FeBiNJ->mBiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPPhHH,'mNEF,EFBiNJ->mBiJ');

Vt3B.oPhH = d1 + d2;
end