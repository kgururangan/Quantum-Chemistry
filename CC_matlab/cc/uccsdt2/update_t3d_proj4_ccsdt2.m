function [t3d] = update_t3d_proj4_ccsdt2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift)

    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

    H1B = HBar_t.H1B; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    [Vt3C] = get_t3d_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys);

    I2C_vvov = H2C.vvov+einsum_kg(H1B.ov,t2c,'me,abim->abie');
    
    % MM23D
    M23_D1 = -einsum_kg(H2C.vooo(PB,:,hB,HB) + Vt3C.PohH,t2c(PB,pB,:,HB),'AmiJ,BcmK->ABciJK'); 
    M23_D1 = M23_D1 - permute(M23_D1,[1,2,3,4,6,5]) - permute(M23_D1,[2,1,3,4,5,6]) + permute(M23_D1,[2,1,3,4,6,5]);
    M23_D2 = einsum_kg(H2C.vooo(pB,:,hB,HB) + Vt3C.pohH,t2c(PB,PB,:,HB),'cmiJ,BAmK->ABciJK');
    M23_D2 = M23_D2 - permute(M23_D2,[1,2,3,4,6,5]);
    M23_D3 = einsum_kg(H2C.vooo(PB,:,HB,HB) + Vt3C.PoHH,t2c(PB,pB,:,hB),'AmKJ,Bcmi->ABciJK');
    M23_D3 = M23_D3 - permute(M23_D3,[2,1,3,4,5,6]);
    M23_D4 = -einsum_kg(H2C.vooo(pB,:,HB,HB) + Vt3C.poHH,t2c(PB,PB,:,hB),'cmKJ,BAmi->ABciJK');
    
    M23_D1234 = M23_D1 + M23_D2 + M23_D3 + M23_D4;
         
    M23_D5 = +einsum_kg(I2C_vvov(PB,PB,hB,:) + Vt3C.PPhv,t2c(:,pB,HB,HB),'ABie,ecJK->ABciJK');
    M23_D6 = -einsum_kg(I2C_vvov(PB,PB,HB,:) + Vt3C.PPHv,t2c(:,pB,hB,HB),'ABJe,eciK->ABciJK');
    M23_D6 = M23_D6 - permute(M23_D6,[1,2,3,4,6,5]);
    M23_D7 = -einsum_kg(I2C_vvov(PB,pB,hB,:) + Vt3C.Pphv,t2c(:,PB,HB,HB),'Acie,eBJK->ABciJK');
    M23_D7 = M23_D7 - permute(M23_D7,[2,1,3,4,5,6]);
    M23_D8 = +einsum_kg(I2C_vvov(PB,pB,HB,:) + Vt3C.PpHv,t2c(:,PB,hB,HB),'AcJe,eBiK->ABciJK');
    M23_D8 = M23_D8 - permute(M23_D8,[2,1,3,4,5,6]) - permute(M23_D8,[1,2,3,4,6,5]) + permute(M23_D8,[2,1,3,4,6,5]);
    
    M23_D5678 = M23_D5 + M23_D6 + M23_D7 + M23_D8;
     
    MM23D = M23_D1234 + M23_D5678;

    % (HBar*T3)_C   
D1 = +einsum_kg(H1B.oo(hB,HB),T3D.PPphhH,'mK,ABcimJ->ABciJK')...
-einsum_kg(H1B.oo(HB,HB),T3D.PPphHH,'MK,ABciJM->ABciJK');
D1 = D1 - permute(D1,[1,2,3,4,6,5]);

D2 = +einsum_kg(H1B.oo(hB,hB),T3D.PPphHH,'mi,ABcmKJ->ABciJK')...
+einsum_kg(H1B.oo(HB,hB),T3D.PPpHHH,'Mi,ABcKJM->ABciJK');

D3 = +einsum_kg(H1B.vv(pB,pB),T3D.PPphHH,'ce,ABeiJK->ABciJK')...
+einsum_kg(H1B.vv(pB,PB),T3D.PPPhHH,'cE,ABEiJK->ABciJK');

D4 = -einsum_kg(H1B.vv(PB,pB),T3D.PpphHH,'Be,AceiJK->ABciJK')...
+einsum_kg(H1B.vv(PB,PB),T3D.PPphHH,'BE,AEciJK->ABciJK');
D4 = D4 - permute(D4,[2,1,3,4,5,6]);

D5 = +0.5*einsum_kg(H2C.oooo(hB,hB,hB,HB),T3D.PPphhH,'mniJ,ABcmnK->ABciJK')...
-einsum_kg(H2C.oooo(HB,hB,hB,HB),T3D.PPphHH,'MniJ,ABcnMK->ABciJK')...
+0.5*einsum_kg(H2C.oooo(HB,HB,hB,HB),T3D.PPpHHH,'MNiJ,ABcMNK->ABciJK');
D5 = D5 - permute(D5,[1,2,3,4,6,5]);

D6 = -0.5*einsum_kg(H2C.oooo(hB,hB,HB,HB),T3D.PPphhh,'mnKJ,ABcmni->ABciJK')...
-einsum_kg(H2C.oooo(HB,hB,HB,HB),T3D.PPphhH,'MnKJ,ABcniM->ABciJK')...
-0.5*einsum_kg(H2C.oooo(HB,HB,HB,HB),T3D.PPphHH,'MNKJ,ABciMN->ABciJK');

D7 = +0.5*einsum_kg(H2C.vvvv(PB,PB,pB,pB),T3D.ppphHH,'ABef,efciJK->ABciJK')...
+einsum_kg(H2C.vvvv(PB,PB,PB,pB),T3D.PpphHH,'ABEf,EfciJK->ABciJK')...
+0.5*einsum_kg(H2C.vvvv(PB,PB,PB,PB),T3D.PPphHH,'ABEF,EFciJK->ABciJK');

D8 = -0.5*einsum_kg(H2C.vvvv(PB,pB,pB,pB),T3D.PpphHH,'Acef,BefiJK->ABciJK')...
+einsum_kg(H2C.vvvv(PB,pB,PB,pB),T3D.PPphHH,'AcEf,EBfiJK->ABciJK')...
-0.5*einsum_kg(H2C.vvvv(PB,pB,PB,PB),T3D.PPPhHH,'AcEF,EFBiJK->ABciJK');
D8 = D8 - permute(D8,[2,1,3,4,5,6]);

D9 = +einsum_kg(H2B.ovvo(hA,PB,pA,hB),T3C.pPphHH,'mAei,eBcmJK->ABciJK')...
+einsum_kg(H2B.ovvo(hA,PB,PA,hB),T3C.PPphHH,'mAEi,EBcmJK->ABciJK')...
+einsum_kg(H2B.ovvo(HA,PB,pA,hB),T3C.pPpHHH,'MAei,eBcMJK->ABciJK')...
+einsum_kg(H2B.ovvo(HA,PB,PA,hB),T3C.PPpHHH,'MAEi,EBcMJK->ABciJK');
D9 = D9 - permute(D9,[2,1,3,4,5,6]);

D10 = -einsum_kg(H2B.ovvo(hA,pB,pA,hB),T3C.pPPhHH,'mcei,eBAmJK->ABciJK')...
-einsum_kg(H2B.ovvo(hA,pB,PA,hB),T3C.PPPhHH,'mcEi,EBAmJK->ABciJK')...
-einsum_kg(H2B.ovvo(HA,pB,pA,hB),T3C.pPPHHH,'Mcei,eBAMJK->ABciJK')...
-einsum_kg(H2B.ovvo(HA,pB,PA,hB),T3C.PPPHHH,'McEi,EBAMJK->ABciJK');

D11 = -einsum_kg(H2B.ovvo(hA,PB,pA,HB),T3C.pPphhH,'mAeJ,eBcmiK->ABciJK')...
-einsum_kg(H2B.ovvo(hA,PB,PA,HB),T3C.PPphhH,'mAEJ,EBcmiK->ABciJK')...
-einsum_kg(H2B.ovvo(HA,PB,pA,HB),T3C.pPpHhH,'MAeJ,eBcMiK->ABciJK')...
-einsum_kg(H2B.ovvo(HA,PB,PA,HB),T3C.PPpHhH,'MAEJ,EBcMiK->ABciJK');
D11 = D11 - permute(D11,[2,1,3,4,5,6]) - permute(D11,[1,2,3,4,6,5]) + permute(D11,[2,1,3,4,6,5]);

D12 = +einsum_kg(H2B.ovvo(hA,pB,pA,HB),T3C.pPPhhH,'mceJ,eBAmiK->ABciJK')...
+einsum_kg(H2B.ovvo(hA,pB,PA,HB),T3C.PPPhhH,'mcEJ,EBAmiK->ABciJK')...
+einsum_kg(H2B.ovvo(HA,pB,pA,HB),T3C.pPPHhH,'MceJ,eBAMiK->ABciJK')...
+einsum_kg(H2B.ovvo(HA,pB,PA,HB),T3C.PPPHhH,'McEJ,EBAMiK->ABciJK');
D12 = D12 - permute(D12,[1,2,3,4,6,5]);

D13 = -einsum_kg(H2C.voov(PB,hB,hB,pB),T3D.PpphHH,'Amie,BecmJK->ABciJK')...
+einsum_kg(H2C.voov(PB,hB,hB,PB),T3D.PPphHH,'AmiE,EBcmJK->ABciJK')...
-einsum_kg(H2C.voov(PB,HB,hB,pB),T3D.PppHHH,'AMie,BecMJK->ABciJK')...
+einsum_kg(H2C.voov(PB,HB,hB,PB),T3D.PPpHHH,'AMiE,EBcMJK->ABciJK');
D13 = D13 - permute(D13,[2,1,3,4,5,6]);

D14 = -einsum_kg(H2C.voov(pB,hB,hB,pB),T3D.PPphHH,'cmie,BAemJK->ABciJK')...
-einsum_kg(H2C.voov(pB,hB,hB,PB),T3D.PPPhHH,'cmiE,EBAmJK->ABciJK')...
-einsum_kg(H2C.voov(pB,HB,hB,pB),T3D.PPpHHH,'cMie,BAeMJK->ABciJK')...
-einsum_kg(H2C.voov(pB,HB,hB,PB),T3D.PPPHHH,'cMiE,EBAMJK->ABciJK');

D15 = +einsum_kg(H2C.voov(PB,hB,HB,pB),T3D.PpphhH,'AmJe,BecmiK->ABciJK')...
-einsum_kg(H2C.voov(PB,hB,HB,PB),T3D.PPphhH,'AmJE,EBcmiK->ABciJK')...
-einsum_kg(H2C.voov(PB,HB,HB,pB),T3D.PpphHH,'AMJe,BeciMK->ABciJK')...
+einsum_kg(H2C.voov(PB,HB,HB,PB),T3D.PPphHH,'AMJE,EBciMK->ABciJK');
D15 = D15 - permute(D15,[2,1,3,4,5,6]) - permute(D15,[1,2,3,4,6,5]) + permute(D15,[2,1,3,4,6,5]);

D16 = +einsum_kg(H2C.voov(pB,hB,HB,pB),T3D.PPphhH,'cmJe,BAemiK->ABciJK')...
+einsum_kg(H2C.voov(pB,hB,HB,PB),T3D.PPPhhH,'cmJE,EBAmiK->ABciJK')...
-einsum_kg(H2C.voov(pB,HB,HB,pB),T3D.PPphHH,'cMJe,BAeiMK->ABciJK')...
-einsum_kg(H2C.voov(pB,HB,HB,PB),T3D.PPPhHH,'cMJE,EBAiMK->ABciJK');
D16 = D16 - permute(D16,[1,2,3,4,6,5]);


    X3D_PPphHH = MM23D + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14 + D15 + D16;

    t3d = T3D.PPphHH;
    for a = 1:sys.Nact_p_beta
        for b = a+1:sys.Nact_p_beta
            for c = 1:sys.Nunact_p_beta
                for i = 1:sys.Nunact_h_beta
                    for j = 1:sys.Nact_h_beta
                        for k = j+1:sys.Nact_h_beta
                            
                            denom = sys.fb_hh(i,i) + sys.fb_HH(j,j) + sys.fb_HH(k,k) - sys.fb_PP(a,a) - sys.fb_PP(b,b) - sys.fb_pp(c,c);
                            
                            t3d(a,b,c,i,j,k) = t3d(a,b,c,i,j,k) + X3D_PPphHH(a,b,c,i,j,k)/(denom-shift);
                            t3d(b,a,c,i,j,k) = -t3d(a,b,c,i,j,k);
                            t3d(a,b,c,i,k,j) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,i,k,j) = t3d(a,b,c,i,j,k);

                        end
                    end
                end
            end
        end
    end
    %T3A.PPphHH = t3a;
       

end

function [Vt3C] = get_t3d_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys)

%     % h1B(mi)
%     H1B.oo = sys.fb_oo ...
%              +einsum_kg(sys.fb_ov,t1b,'me,ei->mi')...
%              +einsum_kg(sys.vB_oovo,t1a,'nmfi,fn->mi')...
%              +einsum_kg(sys.vC_ooov,t1b,'mnif,fn->mi')...
%              +einsum_kg(sys.vB_oovv,t2b,'nmfe,feni->mi')...
%              +0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,efin->mi')...
%              +einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t1b,'me,ei->mi')...
%              +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t1b,'me,ei->mi');
%          
%     % h1B(ae)
%     H1B.vv = sys.fb_vv ...
%              -einsum_kg(sys.fb_ov,t1b,'me,am->ae')...
%              +einsum_kg(sys.vB_ovvv,t1a,'nafe,fn->ae')...
%              +einsum_kg(sys.vC_vovv,t1b,'anef,fn->ae')...
%              -einsum_kg(sys.vB_oovv,t2b,'nmfe,fanm->ae')...
%              -0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,afmn->ae')...
%              -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t1b,'me,am->ae')...
%              -einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t1b,'me,am->ae');
%          
%     % h1B(me)
%     H1B.ov = sys.fb_ov + einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me') + einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
%     
%     % h2B(maei)
%     H2B.ovvo = sys.vB_ovvo ...
%                +einsum_kg(sys.vB_ovvv,t1b,'maef,fi->maei')...
%                -einsum_kg(sys.vB_oovo,t1b,'mnei,an->maei')...
%                -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->maei')...
%                +einsum_kg(sys.vB_oovv,t2c,'mnef,afin->maei')...
%                +einsum_kg(sys.vA_oovv,t2b,'mnef,fani->maei');
% 
%     % h2C(amij)
%     H2C.vooo =  sys.vC_vooo ...
% 	    	+einsum_kg(sys.vC_vovo,t1b,'amej,ei->amij')...
% 	    			-einsum_kg(sys.vC_vovo,t1b,'amei,ej->amij')...
% 	    	-einsum_kg(sys.vC_oooo,t1b,'nmij,an->amij')...
% 	    	+einsum_kg(sys.vB_oovo,t2b,'nmfj,fani->amij')...
% 	    			-einsum_kg(sys.vB_oovo,t2b,'nmfi,fanj->amij')...
% 	    	+0.5*einsum_kg(sys.vC_vovv,t2c,'amfe,feij->amij')...
% 	    	+einsum_kg(sys.vC_ooov,t2c,'mnjf,afin->amij')...
% 	    			-einsum_kg(sys.vC_ooov,t2c,'mnif,afjn->amij')...
% 	    	+einsum_kg(sys.fb_ov,t2c,'me,aeij->amij')...
% 	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t2c,'me,aeij->amij')...
% 	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,ej->nmfj'),t2b,'nmfj,fani->amij')...
% 	    			-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,ei->nmfi'),t2b,'nmfi,fanj->amij')...
% 	    	-einsum_kg(einsum_kg(sys.vC_ooov,t1b,'mnjf,an->majf'),t1b,'majf,fi->amij')...
% 	    			+einsum_kg(einsum_kg(sys.vC_ooov,t1b,'mnif,an->maif'),t1b,'maif,fj->amij')...
% 	    	+einsum_kg(einsum_kg(sys.vC_vovv,t1b,'amfe,fi->amie'),t1b,'amie,ej->amij')...
% 	    	-0.5*einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,an->maef'),t2c,'maef,feij->amij')...
% 	    	+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ej->mnjf'),t2c,'mnjf,afin->amij')...
% 	    			-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ei->mnif'),t2c,'mnif,afjn->amij')...
% 	    	+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t2c,'me,aeij->amij')...
% 	    	-einsum_kg(einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ej->mnjf'),t1b,'mnjf,fi->mnji'),t1b,'mnji,an->amij');
%     H2C.ovoo = -permute(H2C.vooo,[2,1,3,4]);
% 
%     % h2C(abie)
%     H2C.vvov =  sys.vC_vvov...
% 	    	+einsum_kg(sys.vC_vvvv,t1b,'abfe,fi->abie')...
% 	    	-einsum_kg(sys.vC_ovov,t1b,'nbie,an->abie')...
% 	    			+einsum_kg(sys.vC_ovov,t1b,'naie,bn->abie')...
% 	    	+einsum_kg(sys.vB_ovvv,t2b,'nbfe,fani->abie')...
% 	    			-einsum_kg(sys.vB_ovvv,t2b,'nafe,fbni->abie')...
% 	    	+0.5*einsum_kg(sys.vC_oovo,t2c,'mnei,abnm->abie')...
% 	    	+einsum_kg(sys.vC_ovvv,t2c,'nbfe,afin->abie')...
% 	    			-einsum_kg(sys.vC_ovvv,t2c,'nafe,bfin->abie')...
% 	    	-einsum_kg(sys.fb_ov,t2c,'me,abim->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t2c,'me,abim->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,bm->nbfe'),t2b,'nbfe,fani->abie')...
% 	    			+einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,am->nafe'),t2b,'nafe,fbni->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vC_ovvv,t1b,'nbfe,an->abfe'),t1b,'abfe,fi->abie')...
% 	    			+einsum_kg(einsum_kg(sys.vC_ovvv,t1b,'nafe,bn->bafe'),t1b,'bafe,fi->abie')...
% 	    	+einsum_kg(einsum_kg(sys.vC_oovo,t1b,'mnei,bm->bnei'),t1b,'bnei,an->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bm->bnef'),t2c,'bnef,afin->abie')...
% 	    			+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,am->anef'),t2c,'anef,bfin->abie')...
% 	    	+0.5*einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fi->mnei'),t2c,'mnei,abnm->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t2c,'me,abim->abie')...
% 	    	+einsum_kg(einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bm->bnef'),t1b,'bnef,an->baef'),t1b,'baef,fi->abie');
%      H2C.vvvo = -permute(H2C.vvov,[1,2,4,3]);
%      
%      % h2C(mnij)
%      H2C.oooo = sys.vC_oooo...
% 	      	+einsum_kg(sys.vC_ooov,t1b,'mnie,ej->mnij')...
% 	    			-einsum_kg(sys.vC_ooov,t1b,'mnje,ei->mnij')...
% 	    	+0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,efij->mnij')...
% 	    	+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ei->mnif'),t1b,'mnif,fj->mnij'); 
% 
%      % h2C(abef)
%      H2C.vvvv = sys.vC_vvvv...
% 	     	-einsum_kg(sys.vC_ovvv,t1b,'mbef,am->abef')...
% 	     			+einsum_kg(sys.vC_ovvv,t1b,'maef,bm->abef')...
% 	     	+0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,abmn->abef')...
% 	     	+einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bn->mbef'),t1b,'mbef,am->abef');
% 	      
%      % h2C(amie)
%      H2C.voov = sys.vC_voov ...
%                 -einsum_kg(sys.vC_oovo,t1b,'mnei,an->amie')...
%                 +einsum_kg(sys.vC_vovv,t1b,'amfe,fi->amie')...
%                 -einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->amie')...
%                 +einsum_kg(sys.vC_oovv,t2c,'mnef,afin->amie')...
%                 +einsum_kg(sys.vB_oovv,t2b,'nmfe,fani->amie');
            
    % VT3 intermediates
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

d1 = +einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3C.pPphhH,'nmfe,fAeniJ->AmiJ')...
+einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3C.pPPhhH,'nmfE,fAEniJ->AmiJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3C.PPphhH,'nmFe,FAeniJ->AmiJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3C.PPPhhH,'nmFE,FAEniJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3C.pPpHhH,'Nmfe,fAeNiJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3C.pPPHhH,'NmfE,fAENiJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3C.PPpHhH,'NmFe,FAeNiJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPPHhH,'NmFE,FAENiJ->AmiJ');

d2 = -0.5*einsum_kg(sys.vC_oovv(:,hB,pB,pB),T3D.PpphhH,'mnef,AefinJ->AmiJ')...
+einsum_kg(sys.vC_oovv(:,hB,pB,PB),T3D.PPphhH,'mneF,AFeinJ->AmiJ')...
-0.5*einsum_kg(sys.vC_oovv(:,hB,PB,PB),T3D.PPPhhH,'mnEF,AEFinJ->AmiJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,pB,pB),T3D.PpphHH,'mNef,AefiJN->AmiJ')...
-einsum_kg(sys.vC_oovv(:,HB,pB,PB),T3D.PPphHH,'mNeF,AFeiJN->AmiJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPPhHH,'mNEF,AEFiJN->AmiJ');

Vt3C.PohH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3C.ppphhH,'nmfe,faeniJ->amiJ')...
-einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3C.pPphhH,'nmfE,fEaniJ->amiJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3C.PpphhH,'nmFe,FaeniJ->amiJ')...
-einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3C.PPphhH,'nmFE,FEaniJ->amiJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3C.pppHhH,'Nmfe,faeNiJ->amiJ')...
-einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3C.pPpHhH,'NmfE,fEaNiJ->amiJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3C.PppHhH,'NmFe,FaeNiJ->amiJ')...
-einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPpHhH,'NmFE,FEaNiJ->amiJ');

d2 = -0.5*einsum_kg(sys.vC_oovv(:,hB,pB,pB),T3D.ppphhH,'mnef,aefinJ->amiJ')...
-einsum_kg(sys.vC_oovv(:,hB,pB,PB),T3D.PpphhH,'mneF,FaeinJ->amiJ')...
-0.5*einsum_kg(sys.vC_oovv(:,hB,PB,PB),T3D.PPphhH,'mnEF,EFainJ->amiJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,pB,pB),T3D.ppphHH,'mNef,aefiJN->amiJ')...
+einsum_kg(sys.vC_oovv(:,HB,pB,PB),T3D.PpphHH,'mNeF,FaeiJN->amiJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPphHH,'mNEF,EFaiJN->amiJ');

Vt3C.pohH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3C.pPphHH,'nmfe,fAenIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3C.pPPhHH,'nmfE,fAEnIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3C.PPphHH,'nmFe,FAenIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3C.PPPhHH,'nmFE,FAEnIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3C.pPpHHH,'Nmfe,fAeNIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3C.pPPHHH,'NmfE,fAENIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3C.PPpHHH,'NmFe,FAeNIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPPHHH,'NmFE,FAENIJ->AmIJ');

d2 = +0.5*einsum_kg(sys.vC_oovv(:,hB,pB,pB),T3D.PpphHH,'mnef,AefnIJ->AmIJ')...
-einsum_kg(sys.vC_oovv(:,hB,pB,PB),T3D.PPphHH,'mneF,AFenIJ->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,hB,PB,PB),T3D.PPPhHH,'mnEF,AEFnIJ->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,pB,pB),T3D.PppHHH,'mNef,AefIJN->AmIJ')...
-einsum_kg(sys.vC_oovv(:,HB,pB,PB),T3D.PPpHHH,'mNeF,AFeIJN->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPPHHH,'mNEF,AEFIJN->AmIJ');

Vt3C.PoHH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3C.ppphHH,'nmfe,faenIJ->amIJ')...
-einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3C.pPphHH,'nmfE,fEanIJ->amIJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3C.PpphHH,'nmFe,FaenIJ->amIJ')...
-einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3C.PPphHH,'nmFE,FEanIJ->amIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3C.pppHHH,'Nmfe,faeNIJ->amIJ')...
-einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3C.pPpHHH,'NmfE,fEaNIJ->amIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3C.PppHHH,'NmFe,FaeNIJ->amIJ')...
-einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPpHHH,'NmFE,FEaNIJ->amIJ');

d2 = +0.5*einsum_kg(sys.vC_oovv(:,hB,pB,pB),T3D.ppphHH,'mnef,aefnIJ->amIJ')...
+einsum_kg(sys.vC_oovv(:,hB,pB,PB),T3D.PpphHH,'mneF,FaenIJ->amIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,hB,PB,PB),T3D.PPphHH,'mnEF,EFanIJ->amIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,pB,pB),T3D.pppHHH,'mNef,aefIJN->amIJ')...
+einsum_kg(sys.vC_oovv(:,HB,pB,PB),T3D.PppHHH,'mNeF,FaeIJN->amIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPpHHH,'mNEF,EFaIJN->amIJ');

Vt3C.poHH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3C.pPPhhh,'nmfe,fABnim->ABie')...
+einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3C.PPPhhh,'nmFe,FABnim->ABie')...
+einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3C.pPPHhh,'Nmfe,fABNim->ABie')...
+einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3C.PPPHhh,'NmFe,FABNim->ABie')...
+einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3C.pPPhhH,'nMfe,fABniM->ABie')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3C.PPPhhH,'nMFe,FABniM->ABie')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3C.pPPHhH,'NMfe,fABNiM->ABie')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPPHhH,'NMFe,FABNiM->ABie');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,hB,:,pB),T3D.PPphhh,'mnef,ABfimn->ABie')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,:,PB),T3D.PPPhhh,'mneF,ABFimn->ABie')...
-einsum_kg(sys.vC_oovv(HB,hB,:,pB),T3D.PPphhH,'Mnef,ABfinM->ABie')...
-einsum_kg(sys.vC_oovv(HB,hB,:,PB),T3D.PPPhhH,'MneF,ABFinM->ABie')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,pB),T3D.PPphHH,'MNef,ABfiMN->ABie')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPPhHH,'MNeF,ABFiMN->ABie');

Vt3C.PPhv = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3C.pPPhhH,'nmfe,fABnmI->ABIe')...
-einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3C.PPPhhH,'nmFe,FABnmI->ABIe')...
-einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3C.pPPHhH,'Nmfe,fABNmI->ABIe')...
-einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3C.PPPHhH,'NmFe,FABNmI->ABIe')...
+einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3C.pPPhHH,'nMfe,fABnIM->ABIe')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3C.PPPhHH,'nMFe,FABnIM->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3C.pPPHHH,'NMfe,fABNIM->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPPHHH,'NMFe,FABNIM->ABIe');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,hB,:,pB),T3D.PPphhH,'mnef,ABfmnI->ABIe')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,:,PB),T3D.PPPhhH,'mneF,ABFmnI->ABIe')...
+einsum_kg(sys.vC_oovv(HB,hB,:,pB),T3D.PPphHH,'Mnef,ABfnIM->ABIe')...
+einsum_kg(sys.vC_oovv(HB,hB,:,PB),T3D.PPPhHH,'MneF,ABFnIM->ABIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,pB),T3D.PPpHHH,'MNef,ABfIMN->ABIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPPHHH,'MNeF,ABFIMN->ABIe');

Vt3C.PPHv = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3C.pPphhh,'nmfe,fAbnim->Abie')...
+einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3C.PPphhh,'nmFe,FAbnim->Abie')...
+einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3C.pPpHhh,'Nmfe,fAbNim->Abie')...
+einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3C.PPpHhh,'NmFe,FAbNim->Abie')...
+einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3C.pPphhH,'nMfe,fAbniM->Abie')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3C.PPphhH,'nMFe,FAbniM->Abie')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3C.pPpHhH,'NMfe,fAbNiM->Abie')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPpHhH,'NMFe,FAbNiM->Abie');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,hB,:,pB),T3D.Ppphhh,'mnef,Abfimn->Abie')...
-0.5*einsum_kg(sys.vC_oovv(hB,hB,:,PB),T3D.PPphhh,'mneF,AFbimn->Abie')...
-einsum_kg(sys.vC_oovv(HB,hB,:,pB),T3D.PpphhH,'Mnef,AbfinM->Abie')...
+einsum_kg(sys.vC_oovv(HB,hB,:,PB),T3D.PPphhH,'MneF,AFbinM->Abie')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,pB),T3D.PpphHH,'MNef,AbfiMN->Abie')...
-0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPphHH,'MNeF,AFbiMN->Abie');

Vt3C.Pphv = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3C.pPphhH,'nmfe,fAbnmI->AbIe')...
-einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3C.PPphhH,'nmFe,FAbnmI->AbIe')...
-einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3C.pPpHhH,'Nmfe,fAbNmI->AbIe')...
-einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3C.PPpHhH,'NmFe,FAbNmI->AbIe')...
+einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3C.pPphHH,'nMfe,fAbnIM->AbIe')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3C.PPphHH,'nMFe,FAbnIM->AbIe')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3C.pPpHHH,'NMfe,fAbNIM->AbIe')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPpHHH,'NMFe,FAbNIM->AbIe');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,hB,:,pB),T3D.PpphhH,'mnef,AbfmnI->AbIe')...
-0.5*einsum_kg(sys.vC_oovv(hB,hB,:,PB),T3D.PPphhH,'mneF,AFbmnI->AbIe')...
+einsum_kg(sys.vC_oovv(HB,hB,:,pB),T3D.PpphHH,'Mnef,AbfnIM->AbIe')...
-einsum_kg(sys.vC_oovv(HB,hB,:,PB),T3D.PPphHH,'MneF,AFbnIM->AbIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,pB),T3D.PppHHH,'MNef,AbfIMN->AbIe')...
-0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPpHHH,'MNeF,AFbIMN->AbIe');

Vt3C.PpHv = -d1 - d2;
            
end