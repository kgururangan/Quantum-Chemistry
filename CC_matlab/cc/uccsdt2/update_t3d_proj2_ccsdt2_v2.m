function [t3d] = update_t3d_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift)

    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

    H1B = HBar_t.H1B; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    [Vt3C] = get_t3d_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys);

    I2C_vvov = H2C.vvov+einsum_kg(H1B.ov,t2c,'me,abim->abie');
    
    % MM23D
    M23_D1 =   -einsum_kg(H2C.vooo(PB,:,HB,HB) + Vt3C.PoHH,t2c(PB,pB,:,HB),'AmIJ,BcmK->ABcIJK'); % (ab)(k/ij)
    M23_D1 = M23_D1 - permute(M23_D1,[2,1,3,4,5,6]);

    M23_D2 =   +einsum_kg(H2C.vooo(pB,:,HB,HB) + Vt3C.poHH,t2c(PB,PB,:,HB),'cmIJ,BAmK->ABcIJK'); % (k/ij)
         
    M23_D3 =   +einsum_kg(I2C_vvov(PB,PB,HB,:) + Vt3C.PPHv,t2c(:,pB,HB,HB),'ABIe,ecJK->ABcIJK'); % (i/jk)

    M23_D4 =   -einsum_kg(I2C_vvov(PB,pB,HB,:) + Vt3C.PpHv,t2c(:,PB,HB,HB),'AcIe,eBJK->ABcIJK'); % (ab)(i/jk)
    M23_D4 = M23_D4 - permute(M23_D4,[2,1,3,4,5,6]);

    M23_12 = M23_D1 + M23_D2;
    M23_12 = M23_12 - permute(M23_12,[1,2,3,6,5,4]) - permute(M23_12,[1,2,3,4,6,5]);
    
    M23_34 = M23_D3 + M23_D4;
    M23_34 = M23_34 - permute(M23_34,[1,2,3,5,4,6]) - permute(M23_34,[1,2,3,6,5,4]);

    MM23D = M23_12 + M23_34;

    % (HBar*T3)_C   
D1 = -einsum_kg(H1B.oo(hB,HB),T3D.PPphHH,'mK,ABcmIJ->ABcIJK')...
-einsum_kg(H1B.oo(HB,HB),T3D.PPpHHH,'MK,ABcIJM->ABcIJK');
D1 = D1 - permute(D1,[1,2,3,6,5,4]) - permute(D1,[1,2,3,4,6,5]);

D2 = +einsum_kg(H1B.vv(pB,pB),T3D.PPpHHH,'ce,ABeIJK->ABcIJK')...
+einsum_kg(H1B.vv(pB,PB),T3D.PPPHHH,'cE,ABEIJK->ABcIJK');

D3 = +einsum_kg(H1B.vv(PB,PB),T3D.PPpHHH,'BE,AEcIJK->ABcIJK');
D3 = D3 - permute(D3,[2,1,3,4,5,6]);

D4 = -einsum_kg(H2C.oooo(HB,hB,HB,HB),T3D.PPphHH,'MnIJ,ABcnMK->ABcIJK')...
+0.5*einsum_kg(H2C.oooo(HB,HB,HB,HB),T3D.PPpHHH,'MNIJ,ABcMNK->ABcIJK');
D4 = D4 - permute(D4,[1,2,3,6,5,4]) - permute(D4,[1,2,3,4,6,5]);

D5 = +0.5*einsum_kg(H2C.vvvv(PB,PB,PB,PB),T3D.PPpHHH,'ABEF,EFcIJK->ABcIJK');

D6 = +einsum_kg(H2C.vvvv(PB,pB,PB,pB),T3D.PPpHHH,'AcEf,EBfIJK->ABcIJK')...
-0.5*einsum_kg(H2C.vvvv(PB,pB,PB,PB),T3D.PPPHHH,'AcEF,EFBIJK->ABcIJK');
D6 = D6 - permute(D6,[2,1,3,4,5,6]);

D7 = +einsum_kg(H2B.ovvo(hA,PB,PA,HB),T3C.PPphHH,'mAEI,EBcmJK->ABcIJK')...
+einsum_kg(H2B.ovvo(HA,PB,PA,HB),T3C.PPpHHH,'MAEI,EBcMJK->ABcIJK');
D7 = D7 - permute(D7,[2,1,3,4,5,6]);
D7 = D7 - permute(D7,[1,2,3,5,4,6]) - permute(D7,[1,2,3,6,5,4]);

D8 = -einsum_kg(H2B.ovvo(hA,pB,pA,HB),T3C.pPPhHH,'mceI,eBAmJK->ABcIJK')...
-einsum_kg(H2B.ovvo(hA,pB,PA,HB),T3C.PPPhHH,'mcEI,EBAmJK->ABcIJK')...
-einsum_kg(H2B.ovvo(HA,pB,pA,HB),T3C.pPPHHH,'MceI,eBAMJK->ABcIJK')...
-einsum_kg(H2B.ovvo(HA,pB,PA,HB),T3C.PPPHHH,'McEI,EBAMJK->ABcIJK');
D8 = D8 - permute(D8,[1,2,3,5,4,6]) - permute(D8,[1,2,3,6,5,4]);

D9 = +einsum_kg(H2C.voov(PB,hB,HB,PB),T3D.PPphHH,'AmIE,EBcmJK->ABcIJK')...
+einsum_kg(H2C.voov(PB,HB,HB,PB),T3D.PPpHHH,'AMIE,EBcMJK->ABcIJK');
D9 = D9 - permute(D9,[2,1,3,4,5,6]);
D9 = D9 - permute(D9,[1,2,3,5,4,6]) - permute(D9,[1,2,3,6,5,4]);

D10 = -einsum_kg(H2C.voov(pB,hB,HB,pB),T3D.PPphHH,'cmIe,BAemJK->ABcIJK')...
-einsum_kg(H2C.voov(pB,hB,HB,PB),T3D.PPPhHH,'cmIE,EBAmJK->ABcIJK')...
-einsum_kg(H2C.voov(pB,HB,HB,pB),T3D.PPpHHH,'cMIe,BAeMJK->ABcIJK')...
-einsum_kg(H2C.voov(pB,HB,HB,PB),T3D.PPPHHH,'cMIE,EBAMJK->ABcIJK');
D10 = D10 - permute(D10,[1,2,3,5,4,6]) - permute(D10,[1,2,3,6,5,4]);

    X3D_PPpHHH = MM23D + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10;

%
    t3d = T3D.PPpHHH;
    for a = 1:sys.Nact_p_beta
        for b = a+1:sys.Nact_p_beta
            for c = 1:sys.Nunact_p_beta
                for i = 1:sys.Nact_h_beta
                    for j = i+1:sys.Nact_h_beta
                        for k =j+1:sys.Nact_h_beta
                            denom = sys.fb_HH(i,i) + sys.fb_HH(j,j) + sys.fb_HH(k,k) - sys.fb_PP(a,a) - sys.fb_PP(b,b) - sys.fa_pp(c,c);
                            
                            t3d(a,b,c,i,j,k) = t3d(a,b,c,i,j,k) + X3D_PPpHHH(a,b,c,i,j,k)/(denom-shift);     
                            t3d(a,b,c,k,i,j) = t3d(a,b,c,i,j,k);
                            t3d(a,b,c,j,k,i) = t3d(a,b,c,i,j,k);
                            t3d(a,b,c,i,k,j) = -t3d(a,b,c,i,j,k);
                            t3d(a,b,c,j,i,k) = -t3d(a,b,c,i,j,k);
                            t3d(a,b,c,k,j,i) = -t3d(a,b,c,i,j,k);
                            
                            t3d(b,a,c,i,j,k) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,k,i,j) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,j,k,i) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,i,k,j) = t3d(a,b,c,i,j,k);
                            t3d(b,a,c,j,i,k) = t3d(a,b,c,i,j,k);
                            t3d(b,a,c,k,j,i) = t3d(a,b,c,i,j,k);
                            
                        end
                    end
                end
            end
        end
    end
    %T3A.PPpHHH = t3a;
       

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

d1 = +einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3C.pPPhHH,'nmfE,fAEnIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3C.PPphHH,'nmFe,FAenIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3C.PPPhHH,'nmFE,FAEnIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3C.pPPHHH,'NmfE,fAENIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3C.PPpHHH,'NmFe,FAeNIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPPHHH,'NmFE,FAENIJ->AmIJ');

d2 = -einsum_kg(sys.vC_oovv(:,hB,pB,PB),T3D.PPphHH,'mneF,AFenIJ->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,hB,PB,PB),T3D.PPPhHH,'mnEF,AEFnIJ->AmIJ')...
-einsum_kg(sys.vC_oovv(:,HB,pB,PB),T3D.PPpHHH,'mNeF,AFeIJN->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPPHHH,'mNEF,AEFIJN->AmIJ');

Vt3C.PoHH = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3C.PPphHH,'nmFE,FEanIJ->amIJ')...
-einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPpHHH,'NmFE,FEaNIJ->amIJ');

d2 = +0.5*einsum_kg(sys.vC_oovv(:,hB,PB,PB),T3D.PPphHH,'mnEF,EFanIJ->amIJ')...
+0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPpHHH,'mNEF,EFaIJN->amIJ');

Vt3C.poHH = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3C.pPPHhH,'Nmfe,fABNmI->ABIe')...
-einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3C.PPPHhH,'NmFe,FABNmI->ABIe')...
+einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3C.pPPhHH,'nMfe,fABnIM->ABIe')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3C.PPPhHH,'nMFe,FABnIM->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3C.pPPHHH,'NMfe,fABNIM->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPPHHH,'NMFe,FABNIM->ABIe');

d2 = +einsum_kg(sys.vC_oovv(HB,hB,:,pB),T3D.PPphHH,'Mnef,ABfnIM->ABIe')...
+einsum_kg(sys.vC_oovv(HB,hB,:,PB),T3D.PPPhHH,'MneF,ABFnIM->ABIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,pB),T3D.PPpHHH,'MNef,ABfIMN->ABIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPPHHH,'MNeF,ABFIMN->ABIe');

Vt3C.PPHv = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3C.PPpHhH,'NmFe,FAbNmI->AbIe')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3C.PPphHH,'nMFe,FAbnIM->AbIe')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPpHHH,'NMFe,FAbNIM->AbIe');

d2 = -einsum_kg(sys.vC_oovv(HB,hB,:,PB),T3D.PPphHH,'MneF,AFbnIM->AbIe')...
-0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPpHHH,'MNeF,AFbIMN->AbIe');

Vt3C.PpHv = -d1 - d2;
            
end