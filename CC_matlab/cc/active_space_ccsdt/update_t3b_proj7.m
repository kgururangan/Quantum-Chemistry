function [X3B_PPpHHh] = update_t3b_proj7(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys,shift)

    % store active T3 as structs T3A,B,C,D,...
    % e.g. T3B.PPPHHH, etc.
    
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
    
    [H1A,H1B,H2A,H2B,H2C,Vt3A,Vt3B] = get_t3b_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys);
    
    % MM23B + (V*T2*T3)_C    
    I2B_ovoo = H2B.ovoo - einsum_kg(H1A.ov,t2b,'me,ecjk->mcjk');
    I2B_vooo = H2B.vooo - einsum_kg(H1B.ov,t2b,'me,aeik->amik');
    I2A_vooo = H2A.vooo - einsum_kg(H1A.ov,t2a,'me,aeij->amij');   
    
    M23_D1 = einsum_kg(H2B.vvvo(PA,pB,:,hB) + Vt3B.Ppvh,t2a(PA,:,HA,HA),'Bcek,AeIJ->ABcIJk'); % (ab)
    M23_D2 = -einsum_kg(I2B_ovoo(:,pB,HA,hB) + Vt3B.opHh,t2a(PA,PA,HA,:),'mcJk,ABIm->ABcIJk');
    M23_D3 = +einsum_kg(H2B.vvov(PA,pB,HA,:) + Vt3B.PpHv,t2b(PA,:,HA,hB),'AcIe,BeJk->ABcIJk');
    M23_D4 = -einsum_kg(I2B_vooo(PA,:,HA,hB) + Vt3B.PoHh,t2b(PA,pB,HA,:),'AmIk,BcJm->ABcIJk');
    M23_D5 = +einsum_kg(H2A.vvov(PA,PA,HA,:) + Vt3A.PPHv,t2b(:,pB,HA,hB),'ABIe,ecJk->ABcIJk');
    M23_D6 = -einsum_kg(I2A_vooo(PA,:,HA,HA) + Vt3A.PoHH,t2b(PA,pB,:,hB),'AmIJ,Bcmk->ABcIJk');
    
    M23_34 = M23_D3 + M23_D4;
    M23_34 = M23_34 - permute(M23_34,[2,1,3,4,5,6]) - permute(M23_34,[1,2,3,5,4,6]) + permute(M23_34,[2,1,3,5,4,6]);
    M23_25 = M23_D2 + M23_D5;
    M23_25 = M23_25 - permute(M23_25,[1,2,3,5,4,6]);
    M23_16 = M23_D1 + M23_D6;
    M23_16 = M23_16 - permute(M23_16,[2,1,3,4,5,6]);
        
    MM23B = M23_16 + M23_25 + M23_34;
    
 
    %%%%%

    
D1 = -einsum_kg(H1A.oo(hA,HA),T3B.PPphHh,'mI,ABcmJk->ABcIJk')...
-einsum_kg(H1A.oo(HA,HA),T3B.PPpHHh,'MI,ABcMJk->ABcIJk');
D1 = D1 - permute(D1,[1,2,3,5,4,6]);

D2 = -einsum_kg(H1B.oo(hB,hB),T3B.PPpHHh,'mk,ABcIJm->ABcIJk')...
-einsum_kg(H1B.oo(HB,hB),T3B.PPpHHH,'Mk,ABcIJM->ABcIJk');

D3 = -einsum_kg(H1A.vv(PA,pA),T3B.PppHHh,'Ae,BecIJk->ABcIJk')...
+einsum_kg(H1A.vv(PA,PA),T3B.PPpHHh,'AE,EBcIJk->ABcIJk');
D3 = D3 - permute(D3,[2,1,3,4,5,6]);

D4 = +einsum_kg(H1B.vv(pB,pB),T3B.PPpHHh,'ce,ABeIJk->ABcIJk')...
+einsum_kg(H1B.vv(pB,PB),T3B.PPPHHh,'cE,ABEIJk->ABcIJk');

D5 = +0.5*einsum_kg(H2A.oooo(hA,hA,HA,HA),T3B.PPphhh,'mnIJ,ABcmnk->ABcIJk')...
-einsum_kg(H2A.oooo(HA,hA,HA,HA),T3B.PPphHh,'MnIJ,ABcnMk->ABcIJk')...
+0.5*einsum_kg(H2A.oooo(HA,HA,HA,HA),T3B.PPpHHh,'MNIJ,ABcMNk->ABcIJk');

D6 = -einsum_kg(H2B.oooo(hA,hB,HA,hB),T3B.PPphHh,'mnJk,ABcmIn->ABcIJk')...
-einsum_kg(H2B.oooo(hA,HB,HA,hB),T3B.PPphHH,'mNJk,ABcmIN->ABcIJk')...
+einsum_kg(H2B.oooo(HA,hB,HA,hB),T3B.PPpHHh,'MnJk,ABcIMn->ABcIJk')...
+einsum_kg(H2B.oooo(HA,HB,HA,hB),T3B.PPpHHH,'MNJk,ABcIMN->ABcIJk');
D6 = D6 - permute(D6,[1,2,3,5,4,6]);

D7 = +0.5*einsum_kg(H2A.vvvv(PA,PA,pA,pA),T3B.pppHHh,'ABef,efcIJk->ABcIJk')...
+einsum_kg(H2A.vvvv(PA,PA,PA,pA),T3B.PppHHh,'ABEf,EfcIJk->ABcIJk')...
+0.5*einsum_kg(H2A.vvvv(PA,PA,PA,PA),T3B.PPpHHh,'ABEF,EFcIJk->ABcIJk');

D8 = +einsum_kg(H2B.vvvv(PA,pB,pA,pB),T3B.PppHHh,'Bcef,AefIJk->ABcIJk')...
+einsum_kg(H2B.vvvv(PA,pB,pA,PB),T3B.PpPHHh,'BceF,AeFIJk->ABcIJk')...
+einsum_kg(H2B.vvvv(PA,pB,PA,pB),T3B.PPpHHh,'BcEf,AEfIJk->ABcIJk')...
+einsum_kg(H2B.vvvv(PA,pB,PA,PB),T3B.PPPHHh,'BcEF,AEFIJk->ABcIJk');
D8 = D8 - permute(D8,[2,1,3,4,5,6]);

D9 = -einsum_kg(H2A.voov(PA,hA,HA,pA),T3B.PpphHh,'AmIe,BecmJk->ABcIJk')...
+einsum_kg(H2A.voov(PA,hA,HA,PA),T3B.PPphHh,'AmIE,EBcmJk->ABcIJk')...
-einsum_kg(H2A.voov(PA,HA,HA,pA),T3B.PppHHh,'AMIe,BecMJk->ABcIJk')...
+einsum_kg(H2A.voov(PA,HA,HA,PA),T3B.PPpHHh,'AMIE,EBcMJk->ABcIJk');
D9 = D9 - permute(D9,[2,1,3,4,5,6]);
D9 = D9 - permute(D9,[1,2,3,5,4,6]);

D10 = +einsum_kg(H2B.voov(PA,hB,HA,pB),T3C.PppHhh,'AmIe,BceJkm->ABcIJk')...
-einsum_kg(H2B.voov(PA,hB,HA,PB),T3C.PPpHhh,'AmIE,BEcJkm->ABcIJk')...
+einsum_kg(H2B.voov(PA,HB,HA,pB),T3C.PppHhH,'AMIe,BceJkM->ABcIJk')...
-einsum_kg(H2B.voov(PA,HB,HA,PB),T3C.PPpHhH,'AMIE,BEcJkM->ABcIJk');
D10 = D10 - permute(D10,[2,1,3,4,5,6]);
D10 = D10 - permute(D10,[1,2,3,5,4,6]);

D11 = +einsum_kg(H2B.ovvo(hA,pB,pA,hB),T3A.PPphHH,'mcek,ABemIJ->ABcIJk')...
+einsum_kg(H2B.ovvo(hA,pB,PA,hB),T3A.PPPhHH,'mcEk,ABEmIJ->ABcIJk')...
+einsum_kg(H2B.ovvo(HA,pB,pA,hB),T3A.PPpHHH,'Mcek,ABeIJM->ABcIJk')...
+einsum_kg(H2B.ovvo(HA,pB,PA,hB),T3A.PPPHHH,'McEk,ABEIJM->ABcIJk');

D12 = +einsum_kg(H2C.voov(pB,hB,hB,pB),T3B.PPpHHh,'cmke,ABeIJm->ABcIJk')...
+einsum_kg(H2C.voov(pB,hB,hB,PB),T3B.PPPHHh,'cmkE,ABEIJm->ABcIJk')...
+einsum_kg(H2C.voov(pB,HB,hB,pB),T3B.PPpHHH,'cMke,ABeIJM->ABcIJk')...
+einsum_kg(H2C.voov(pB,HB,hB,PB),T3B.PPPHHH,'cMkE,ABEIJM->ABcIJk');

D13 = +einsum_kg(H2B.vovo(PA,hB,pA,hB),T3B.PppHHh,'Amek,BecIJm->ABcIJk')...
-einsum_kg(H2B.vovo(PA,hB,PA,hB),T3B.PPpHHh,'AmEk,EBcIJm->ABcIJk')...
+einsum_kg(H2B.vovo(PA,HB,pA,hB),T3B.PppHHH,'AMek,BecIJM->ABcIJk')...
-einsum_kg(H2B.vovo(PA,HB,PA,hB),T3B.PPpHHH,'AMEk,EBcIJM->ABcIJk');
D13 = D13 - permute(D13,[2,1,3,4,5,6]);

D14 = -einsum_kg(H2B.ovov(hA,pB,HA,pB),T3B.PPphHh,'mcIe,ABemJk->ABcIJk')...
-einsum_kg(H2B.ovov(hA,pB,HA,PB),T3B.PPPhHh,'mcIE,ABEmJk->ABcIJk')...
-einsum_kg(H2B.ovov(HA,pB,HA,pB),T3B.PPpHHh,'McIe,ABeMJk->ABcIJk')...
-einsum_kg(H2B.ovov(HA,pB,HA,PB),T3B.PPPHHh,'McIE,ABEMJk->ABcIJk');
D14 = D14 - permute(D14,[1,2,3,5,4,6]);


    X3B_PPpHHh = MM23B + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14;

end

function [H1A,H1B,H2A,H2B,H2C,Vt3A,Vt3B,Vt3C] = get_t3b_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys)

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
    
    % h2A(amij)
    H2A.vooo    = sys.vA_vooo ...
                  +einsum_kg(sys.vA_voov,t1a,'amie,ej->amij')...
	    			-einsum_kg(sys.vA_voov,t1a,'amje,ei->amij')...
                  -einsum_kg(sys.vA_oooo,t1a,'nmij,an->amij')...
                  +einsum_kg(sys.fa_ov,t2a,'me,aeij->amij')...
                  +0.5*einsum_kg(sys.vA_vovv,t2a,'amef,efij->amij')...
                  +einsum_kg(sys.vA_oovo,t2a,'nmej,aein->amij')...
	    			-einsum_kg(sys.vA_oovo,t2a,'nmei,aejn->amij')...
                  +einsum_kg(sys.vB_ooov,t2b,'mnje,aein->amij')...
	    			-einsum_kg(sys.vB_ooov,t2b,'mnie,aejn->amij')...
                  -einsum_kg(einsum_kg(sys.vA_oovo,t1a,'nmej,ei->nmij'),t1a,'nmij,an->amij')...
	    			+einsum_kg(einsum_kg(sys.vA_oovo,t1a,'nmei,ej->nmij'),t1a,'nmij,an->amij')...
                  +einsum_kg(einsum_kg(sys.vA_vovv,t1a,'amef,ei->amif'),t1a,'amif,fj->amij')...
                  -0.5*einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmef,an->amef'),t2a,'amef,efij->amij')...
                  +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ej->mnjf'),t2a,'mnjf,afin->amij')...
	    			-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ei->mnif'),t2a,'mnif,afjn->amij')...
                  +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t2a,'me,aeij->amij')...
                  +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ej->mnjf'),t2b,'mnjf,afin->amij')...
	    			-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ei->mnif'),t2b,'mnif,afjn->amij')...
                  +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t2a,'me,aeij->amij')...
                  -einsum_kg(einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmef,fj->nmej'),t1a,'nmej,an->amej'),t1a,'amej,ei->amij');

    % h(abie)
    H2A.vvov  = sys.vA_vvov ...
	    	+einsum_kg(sys.vA_vvvv,t1a,'abfe,fi->abie')...
	    	-einsum_kg(sys.vA_ovov,t1a,'nbie,an->abie')...
	    			+einsum_kg(sys.vA_ovov,t1a,'naie,bn->abie')...
	    	-einsum_kg(sys.fa_ov,t2a,'me,abim->abie')...
	    	+einsum_kg(sys.vA_ovvv,t2a,'nbfe,afin->abie')...
	    			-einsum_kg(sys.vA_ovvv,t2a,'nafe,bfin->abie')...
	    	+0.5*einsum_kg(sys.vA_oovo,t2a,'mnei,abnm->abie')...
	    	+einsum_kg(sys.vB_vovv,t2b,'bnef,afin->abie')...
	    			-einsum_kg(sys.vB_vovv,t2b,'anef,bfin->abie')...
	    	+einsum_kg(einsum_kg(sys.vA_oovo,t1a,'mnei,an->maei'),t1a,'maei,bm->abie')...
	    	-einsum_kg(einsum_kg(sys.vA_vovv,t1a,'bnef,fi->bnei'),t1a,'bnei,an->abie')...
	    			+einsum_kg(einsum_kg(sys.vA_vovv,t1a,'anef,fi->anei'),t1a,'anei,bn->abie')...
	    	+0.5*einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fi->mnei'),t2a,'mnei,abnm->abie')...
	    	-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bm->bnef'),t2a,'bnef,afin->abie')...
	    			+einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,am->anef'),t2a,'anef,bfin->abie')...
	    	-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t2a,'me,abim->abie')...
	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,bm->bnef'),t2b,'bnef,afin->abie')...
	    			+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,am->anef'),t2b,'anef,bfin->abie')...
	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t2a,'me,abim->abie')...
	    	+einsum_kg(einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bm->bnef'),t1a,'bnef,an->baef'),t1a,'baef,fi->abie');

     % h2A(mnij)
     H2A.oooo = sys.vA_oooo...
	     	+einsum_kg(sys.vA_oovo,t1a,'mnej,ei->mnij')...
	     			-einsum_kg(sys.vA_oovo,t1a,'mnei,ej->mnij')...
	     	+0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efij->mnij')...
	     	+einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ei->mnif'),t1a,'mnif,fj->mnij');

     % h2A(abef)
     H2A.vvvv = sys.vA_vvvv...
	     	-einsum_kg(sys.vA_ovvv,t1a,'mbef,am->abef')...
	     			+einsum_kg(sys.vA_ovvv,t1a,'maef,bm->abef')...
	     	+0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,abmn->abef')...
	     	+einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bn->mbef'),t1a,'mbef,am->abef');
        
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
    
d1 = +einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3C.Ppphhh,'mnef,Afbmnj->Abej')...
+einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3C.PPphhh,'mneF,AFbmnj->Abej')...
-einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3C.PpphhH,'mNef,AfbmjN->Abej')...
-einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3C.PPphhH,'mNeF,AFbmjN->Abej')...
+einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3C.PppHhh,'Mnef,AfbMnj->Abej')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3C.PPpHhh,'MneF,AFbMnj->Abej')...
-einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3C.PppHhH,'MNef,AfbMjN->Abej')...
-einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.PPpHhH,'MNeF,AFbMjN->Abej');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3B.Ppphhh,'mnef,Afbmnj->Abej')...
+0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3B.PPphhh,'mneF,AFbmnj->Abej')...
-einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3B.PpphHh,'Mnef,AfbnMj->Abej')...
-einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3B.PPphHh,'MneF,AFbnMj->Abej')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3B.PppHHh,'MNef,AfbMNj->Abej')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PPpHHh,'MNeF,AFbMNj->Abej');

Vt3B.Ppvh = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3C.pppHhh,'mnef,efbInj->mbIj')...
+einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3C.PppHhh,'mnEf,EfbInj->mbIj')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3C.pPpHhh,'mneF,eFbInj->mbIj')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3C.PPpHhh,'mnEF,EFbInj->mbIj')...
-einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3C.pppHhH,'mNef,efbIjN->mbIj')...
-einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PppHhH,'mNEf,EfbIjN->mbIj')...
-einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPpHhH,'mNeF,eFbIjN->mbIj')...
-einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPpHhH,'mNEF,EFbIjN->mbIj');

d2 = -0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3B.ppphHh,'mnef,efbnIj->mbIj')...
+einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3B.PpphHh,'mneF,FebnIj->mbIj')...
-0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3B.PPphHh,'mnEF,EFbnIj->mbIj')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3B.pppHHh,'mNef,efbINj->mbIj')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PppHHh,'mNeF,FebINj->mbIj')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPpHHh,'mNEF,EFbINj->mbIj');

Vt3B.opHh = d1 + d2;
              

d1 = -einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3B.PpphHh,'nmfe,AfbnIm->AbIe')...
-einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3B.PPphHh,'nmFe,AFbnIm->AbIe')...
+einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3B.PppHHh,'Nmfe,AfbINm->AbIe')...
+einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3B.PPpHHh,'NmFe,AFbINm->AbIe')...
-einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3B.PpphHH,'nMfe,AfbnIM->AbIe')...
-einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3B.PPphHH,'nMFe,AFbnIM->AbIe')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3B.PppHHH,'NMfe,AfbINM->AbIe')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3B.PPpHHH,'NMFe,AFbINM->AbIe');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,hB,pB,:),T3C.PppHhh,'nmfe,AfbInm->AbIe')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,PB,:),T3C.PPpHhh,'nmFe,AFbInm->AbIe')...
+einsum_kg(sys.vC_oovv(hB,HB,pB,:),T3C.PppHhH,'nMfe,AfbInM->AbIe')...
+einsum_kg(sys.vC_oovv(hB,HB,PB,:),T3C.PPpHhH,'nMFe,AFbInM->AbIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,pB,:),T3C.PppHHH,'NMfe,AfbINM->AbIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),T3C.PPpHHH,'NMFe,AFbINM->AbIe');

Vt3B.PpHv = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3B.PpphHh,'nmfe,AfenIj->AmIj')...
-einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3B.PpPhHh,'nmfE,AfEnIj->AmIj')...
-einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3B.PPphHh,'nmFe,AFenIj->AmIj')...
-einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3B.PPPhHh,'nmFE,AFEnIj->AmIj')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3B.PppHHh,'Nmfe,AfeINj->AmIj')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3B.PpPHHh,'NmfE,AfEINj->AmIj')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3B.PPpHHh,'NmFe,AFeINj->AmIj')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3B.PPPHHh,'NmFE,AFEINj->AmIj');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,:,pB,pB),T3C.PppHhh,'nmfe,AfeInj->AmIj')...
+einsum_kg(sys.vC_oovv(hB,:,PB,pB),T3C.PPpHhh,'nmFe,AFeInj->AmIj')...
+0.5*einsum_kg(sys.vC_oovv(hB,:,PB,PB),T3C.PPPHhh,'nmFE,AFEInj->AmIj')...
-0.5*einsum_kg(sys.vC_oovv(HB,:,pB,pB),T3C.PppHhH,'Nmfe,AfeIjN->AmIj')...
-einsum_kg(sys.vC_oovv(HB,:,PB,pB),T3C.PPpHhH,'NmFe,AFeIjN->AmIj')...
-0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),T3C.PPPHhH,'NmFE,AFEIjN->AmIj');

Vt3B.PoHh = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3B.PPphHh,'mnef,ABfmIn->ABIe')...
-einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3B.PPPhHh,'mneF,ABFmIn->ABIe')...
-einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3B.PPphHH,'mNef,ABfmIN->ABIe')...
-einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3B.PPPhHH,'mNeF,ABFmIN->ABIe')...
+einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3B.PPpHHh,'Mnef,ABfIMn->ABIe')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3B.PPPHHh,'MneF,ABFIMn->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3B.PPpHHH,'MNef,ABfIMN->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3B.PPPHHH,'MNeF,ABFIMN->ABIe');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3A.PPphhH,'mnef,ABfmnI->ABIe')...
+0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3A.PPPhhH,'mneF,ABFmnI->ABIe')...
+einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3A.PPphHH,'Mnef,ABfnIM->ABIe')...
+einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3A.PPPhHH,'MneF,ABFnIM->ABIe')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3A.PPpHHH,'MNef,ABfIMN->ABIe')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3A.PPPHHH,'MNeF,ABFIMN->ABIe');

Vt3A.PPHv = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3B.PppHHh,'mnef,AefIJn->AmIJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3B.PPpHHh,'mnEf,AEfIJn->AmIJ')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3B.PpPHHh,'mneF,AeFIJn->AmIJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3B.PPPHHh,'mnEF,AEFIJn->AmIJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3B.PppHHH,'mNef,AefIJN->AmIJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3B.PPpHHH,'mNEf,AEfIJN->AmIJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3B.PpPHHH,'mNeF,AeFIJN->AmIJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3B.PPPHHH,'mNEF,AEFIJN->AmIJ');

d2 = +0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3A.PpphHH,'mnef,AefnIJ->AmIJ')...
-einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3A.PPphHH,'mneF,AFenIJ->AmIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3A.PPPhHH,'mnEF,AEFnIJ->AmIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3A.PppHHH,'mNef,AefIJN->AmIJ')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3A.PPpHHH,'mNeF,AFeIJN->AmIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3A.PPPHHH,'mNEF,AEFIJN->AmIJ');

Vt3A.PoHH = d1 + d2;
end