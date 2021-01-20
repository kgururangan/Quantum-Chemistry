function [t3b] = update_t3b_proj5_ccsdt2_v2(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,HBar_t,sys,shift)

    % store active T3 as structs T3A,B,C,D,...
    % e.g. T3B.PPPHHH, etc.
    
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
    
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    [Vt3A,Vt3B] = get_t3b_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys);
    
    % MM23B + (V*T2*T3)_C    
    I2B_ovoo = H2B.ovoo - einsum_kg(H1A.ov,t2b,'me,ecjk->mcjk');
    I2B_vooo = H2B.vooo - einsum_kg(H1B.ov,t2b,'me,aeik->amik');
    I2A_vooo = H2A.vooo - einsum_kg(H1A.ov,t2a,'me,aeij->amij');   
    
    M23_D1 = einsum_kg(H2B.vvvo(PA,PB,:,HB) + Vt3B.PPvH,t2a(PA,:,hA,HB),'BCeK,AeiJ->ABCiJK'); % (ab)
    
    M23_D2 = -einsum_kg(I2B_ovoo(:,PB,HA,HB) + Vt3B.oPHH,t2a(PA,PB,hA,:),'mCJK,ABim->ABCiJK');
    M23_D3 = +einsum_kg(I2B_ovoo(:,PB,hA,HB) + Vt3B.oPhH,t2a(PA,PB,HA,:),'mCiK,ABJm->ABCiJK');
    
    M23_D4 = +einsum_kg(H2B.vvov(PA,PB,hA,:) + Vt3B.PPhv,t2b(PA,:,HA,HB),'ACie,BeJK->ABCiJK'); % (ab)
    M23_D5 = -einsum_kg(H2B.vvov(PA,PB,HA,:) + Vt3B.PPHv,t2b(PA,:,hA,HB),'ACJe,BeiK->ABCiJK'); % (ab)
    
    M23_D6 = -einsum_kg(I2B_vooo(PA,:,hA,HB) + Vt3B.PohH,t2b(PA,PB,HA,:),'AmiK,BCJm->ABCiJK'); % (ab)
    M23_D7 = +einsum_kg(I2B_vooo(PA,:,HA,HB) + Vt3B.PoHH,t2b(PA,PB,hA,:),'AmJK,BCim->ABCiJK'); % (ab)
    
    M23_D8 = +einsum_kg(H2A.vvov(PA,PA,hA,:) + Vt3A.PPhv,t2b(:,PB,HA,HB),'ABie,eCJK->ABCiJK');
    M23_D9 = -einsum_kg(H2A.vvov(PA,PA,HA,:) + Vt3A.PPHv,t2b(:,PB,hA,HB),'ABJe,eCiK->ABCiJK');
    
    M23_D10 = -einsum_kg(I2A_vooo(PA,:,hA,HA) + Vt3A.PohH,t2b(PA,PB,:,HB),'AmiJ,BCmK->ABCiJK'); % (ab)
    
    M23_D1456710 = M23_D1 + M23_D4 + M23_D5 + M23_D6 + M23_D7 + M23_D10;
    M23_D1456710 = M23_D1456710 - permute(M23_D1456710,[2,1,3,4,5,6]);
        
    MM23B = M23_D1456710 + M23_D2 + M23_D3 + M23_D8 + M23_D9;
 
    %%%%%

D1 = -einsum_kg(H1A.oo(hA,hA),T3B.PPPhHH,'mi,ABCmJK->ABCiJK')...
-einsum_kg(H1A.oo(HA,hA),T3B.PPPHHH,'Mi,ABCMJK->ABCiJK')...
-einsum_kg(H1A.oo(HA,HA),T3B.PPPhHH,'MJ,ABCiMK->ABCiJK');

D2 = -einsum_kg(H1B.oo(HB,HB),T3B.PPPhHH,'MK,ABCiJM->ABCiJK');

D3 = -einsum_kg(H1A.vv(PA,pA),T3B.PpPhHH,'Ae,BeCiJK->ABCiJK')...
+einsum_kg(H1A.vv(PA,PA),T3B.PPPhHH,'AE,EBCiJK->ABCiJK');
D3 = D3 - permute(D3,[2,1,3,4,5,6]);

D4 = +einsum_kg(H1B.vv(PB,pB),T3B.PPphHH,'Ce,ABeiJK->ABCiJK')...
+einsum_kg(H1B.vv(PB,PB),T3B.PPPhHH,'CE,ABEiJK->ABCiJK');

D5 = -einsum_kg(H2A.oooo(HA,hA,hA,HA),T3B.PPPhHH,'MniJ,ABCnMK->ABCiJK')...
+0.5*einsum_kg(H2A.oooo(HA,HA,hA,HA),T3B.PPPHHH,'MNiJ,ABCMNK->ABCiJK');

D6 = +einsum_kg(H2B.oooo(HA,HB,HA,HB),T3B.PPPhHH,'MNJK,ABCiMN->ABCiJK')...
+einsum_kg(H2B.oooo(hA,HB,hA,HB),T3B.PPPhHH,'mNiK,ABCmJN->ABCiJK')...
-einsum_kg(H2B.oooo(HA,hB,hA,HB),T3B.PPPHHh,'MniK,ABCJMn->ABCiJK')...
-einsum_kg(H2B.oooo(HA,HB,hA,HB),T3B.PPPHHH,'MNiK,ABCJMN->ABCiJK');

D7 = +einsum_kg(H2A.vvvv(PA,PA,PA,pA),T3B.PpPhHH,'ABEf,EfCiJK->ABCiJK')...
+0.5*einsum_kg(H2A.vvvv(PA,PA,PA,PA),T3B.PPPhHH,'ABEF,EFCiJK->ABCiJK');

D8 = +einsum_kg(H2B.vvvv(PA,PB,pA,PB),T3B.PpPhHH,'BCeF,AeFiJK->ABCiJK')...
+einsum_kg(H2B.vvvv(PA,PB,PA,pB),T3B.PPphHH,'BCEf,AEfiJK->ABCiJK')...
+einsum_kg(H2B.vvvv(PA,PB,PA,PB),T3B.PPPhHH,'BCEF,AEFiJK->ABCiJK');
D8 = D8 - permute(D8,[2,1,3,4,5,6]);

D9 = -einsum_kg(H2A.voov(PA,hA,hA,pA),T3B.PpPhHH,'Amie,BeCmJK->ABCiJK')...
+einsum_kg(H2A.voov(PA,hA,hA,PA),T3B.PPPhHH,'AmiE,EBCmJK->ABCiJK')...
-einsum_kg(H2A.voov(PA,HA,hA,pA),T3B.PpPHHH,'AMie,BeCMJK->ABCiJK')...
+einsum_kg(H2A.voov(PA,HA,hA,PA),T3B.PPPHHH,'AMiE,EBCMJK->ABCiJK')...
-einsum_kg(H2A.voov(PA,HA,HA,pA),T3B.PpPhHH,'AMJe,BeCiMK->ABCiJK')...
+einsum_kg(H2A.voov(PA,HA,HA,PA),T3B.PPPhHH,'AMJE,EBCiMK->ABCiJK');
D9 = D9 - permute(D9,[2,1,3,4,5,6]);

D10 = -einsum_kg(H2B.voov(PA,hB,hA,pB),T3C.PPpHhH,'Amie,BCeJmK->ABCiJK')...
-einsum_kg(H2B.voov(PA,hB,hA,PB),T3C.PPPHhH,'AmiE,BCEJmK->ABCiJK')...
+einsum_kg(H2B.voov(PA,HB,hA,pB),T3C.PPpHHH,'AMie,BCeJKM->ABCiJK')...
+einsum_kg(H2B.voov(PA,HB,hA,PB),T3C.PPPHHH,'AMiE,BCEJKM->ABCiJK')...
-einsum_kg(H2B.voov(PA,HB,HA,pB),T3C.PPphHH,'AMJe,BCeiKM->ABCiJK')...
-einsum_kg(H2B.voov(PA,HB,HA,PB),T3C.PPPhHH,'AMJE,BCEiKM->ABCiJK');
D10 = D10 - permute(D10,[2,1,3,4,5,6]);

D11 = +einsum_kg(H2B.ovvo(HA,PB,pA,HB),T3A.PPphHH,'MCeK,ABeiJM->ABCiJK')...
+einsum_kg(H2B.ovvo(HA,PB,PA,HB),T3A.PPPhHH,'MCEK,ABEiJM->ABCiJK');

D12 = +einsum_kg(H2C.voov(PB,HB,HB,pB),T3B.PPphHH,'CMKe,ABeiJM->ABCiJK')...
+einsum_kg(H2C.voov(PB,HB,HB,PB),T3B.PPPhHH,'CMKE,ABEiJM->ABCiJK');

D13 = +einsum_kg(H2B.vovo(PA,HB,pA,HB),T3B.PpPhHH,'AMeK,BeCiJM->ABCiJK')...
-einsum_kg(H2B.vovo(PA,HB,PA,HB),T3B.PPPhHH,'AMEK,EBCiJM->ABCiJK');
D13 = D13 - permute(D13,[2,1,3,4,5,6]);

D14 = -einsum_kg(H2B.ovov(hA,PB,hA,pB),T3B.PPphHH,'mCie,ABemJK->ABCiJK')...
-einsum_kg(H2B.ovov(hA,PB,hA,PB),T3B.PPPhHH,'mCiE,ABEmJK->ABCiJK')...
-einsum_kg(H2B.ovov(HA,PB,hA,pB),T3B.PPpHHH,'MCie,ABeMJK->ABCiJK')...
-einsum_kg(H2B.ovov(HA,PB,hA,PB),T3B.PPPHHH,'MCiE,ABEMJK->ABCiJK')...
-einsum_kg(H2B.ovov(HA,PB,HA,pB),T3B.PPphHH,'MCJe,ABeiMK->ABCiJK')...
-einsum_kg(H2B.ovov(HA,PB,HA,PB),T3B.PPPhHH,'MCJE,ABEiMK->ABCiJK');

    X3B_PPPhHH = MM23B + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14;

    t3b = T3B.PPPhHH;
    for a = 1:sys.Nact_p_alpha
        for b = a+1:sys.Nact_p_alpha
            for c = 1:sys.Nact_p_beta
                for i = 1:sys.Nunact_h_alpha
                    for j = 1:sys.Nact_h_alpha
                        for k = 1:sys.Nact_h_beta
                            denom=sys.fa_hh(i,i)+sys.fa_HH(j,j)+sys.fb_HH(k,k)-sys.fa_PP(a,a)-sys.fa_PP(b,b)-sys.fb_PP(c,c);
                            t3b(a,b,c,i,j,k) = t3b(a,b,c,i,j,k) + X3B_PPPhHH(a,b,c,i,j,k)/(denom-shift);
                            t3b(b,a,c,i,j,k) = -t3b(a,b,c,i,j,k);
                        end
                    end
                end
            end
        end
    end

end

function [Vt3A,Vt3B] = get_t3b_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys)

%     % h1A(mi)
%     H1A.oo = sys.fa_oo ...
%              +einsum_kg(sys.fa_ov,t1a,'me,ei->mi')...
%              +einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi')...
%              +einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi')...
%              +0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi')...
%              +einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi')...
%              +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t1a,'me,ei->mi')...
%              +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t1a,'me,ei->mi');
%          
%     % h1A(ae)
%     H1A.vv = sys.fa_vv ...
%              -einsum_kg(sys.fa_ov,t1a,'me,am->ae')...
%              +einsum_kg(sys.vA_vovv,t1a,'amef,fm->ae')...
%              +einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae')...
%              -0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,afmn->ae')...
%              -einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ae')...
%              -einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t1a,'me,am->ae')...
%              -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t1a,'me,am->ae');
%          
%     % h1A(me)
%     H1A.ov = sys.fa_ov + einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me') + einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');
% 
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
%     % h2A(amij)
%     H2A.vooo    = sys.vA_vooo ...
%                   +einsum_kg(sys.vA_voov,t1a,'amie,ej->amij')...
% 	    			-einsum_kg(sys.vA_voov,t1a,'amje,ei->amij')...
%                   -einsum_kg(sys.vA_oooo,t1a,'nmij,an->amij')...
%                   +einsum_kg(sys.fa_ov,t2a,'me,aeij->amij')...
%                   +0.5*einsum_kg(sys.vA_vovv,t2a,'amef,efij->amij')...
%                   +einsum_kg(sys.vA_oovo,t2a,'nmej,aein->amij')...
% 	    			-einsum_kg(sys.vA_oovo,t2a,'nmei,aejn->amij')...
%                   +einsum_kg(sys.vB_ooov,t2b,'mnje,aein->amij')...
% 	    			-einsum_kg(sys.vB_ooov,t2b,'mnie,aejn->amij')...
%                   -einsum_kg(einsum_kg(sys.vA_oovo,t1a,'nmej,ei->nmij'),t1a,'nmij,an->amij')...
% 	    			+einsum_kg(einsum_kg(sys.vA_oovo,t1a,'nmei,ej->nmij'),t1a,'nmij,an->amij')...
%                   +einsum_kg(einsum_kg(sys.vA_vovv,t1a,'amef,ei->amif'),t1a,'amif,fj->amij')...
%                   -0.5*einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmef,an->amef'),t2a,'amef,efij->amij')...
%                   +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ej->mnjf'),t2a,'mnjf,afin->amij')...
% 	    			-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ei->mnif'),t2a,'mnif,afjn->amij')...
%                   +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t2a,'me,aeij->amij')...
%                   +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ej->mnjf'),t2b,'mnjf,afin->amij')...
% 	    			-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ei->mnif'),t2b,'mnif,afjn->amij')...
%                   +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t2a,'me,aeij->amij')...
%                   -einsum_kg(einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmef,fj->nmej'),t1a,'nmej,an->amej'),t1a,'amej,ei->amij');
% 
%     % h(abie)
%     H2A.vvov  = sys.vA_vvov ...
% 	    	+einsum_kg(sys.vA_vvvv,t1a,'abfe,fi->abie')...
% 	    	-einsum_kg(sys.vA_ovov,t1a,'nbie,an->abie')...
% 	    			+einsum_kg(sys.vA_ovov,t1a,'naie,bn->abie')...
% 	    	-einsum_kg(sys.fa_ov,t2a,'me,abim->abie')...
% 	    	+einsum_kg(sys.vA_ovvv,t2a,'nbfe,afin->abie')...
% 	    			-einsum_kg(sys.vA_ovvv,t2a,'nafe,bfin->abie')...
% 	    	+0.5*einsum_kg(sys.vA_oovo,t2a,'mnei,abnm->abie')...
% 	    	+einsum_kg(sys.vB_vovv,t2b,'bnef,afin->abie')...
% 	    			-einsum_kg(sys.vB_vovv,t2b,'anef,bfin->abie')...
% 	    	+einsum_kg(einsum_kg(sys.vA_oovo,t1a,'mnei,an->maei'),t1a,'maei,bm->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vA_vovv,t1a,'bnef,fi->bnei'),t1a,'bnei,an->abie')...
% 	    			+einsum_kg(einsum_kg(sys.vA_vovv,t1a,'anef,fi->anei'),t1a,'anei,bn->abie')...
% 	    	+0.5*einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fi->mnei'),t2a,'mnei,abnm->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bm->bnef'),t2a,'bnef,afin->abie')...
% 	    			+einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,am->anef'),t2a,'anef,bfin->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t2a,'me,abim->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,bm->bnef'),t2b,'bnef,afin->abie')...
% 	    			+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,am->anef'),t2b,'anef,bfin->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t2a,'me,abim->abie')...
% 	    	+einsum_kg(einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bm->bnef'),t1a,'bnef,an->baef'),t1a,'baef,fi->abie');
% 
%      % h2A(mnij)
%      H2A.oooo = sys.vA_oooo...
% 	     	+einsum_kg(sys.vA_oovo,t1a,'mnej,ei->mnij')...
% 	     			-einsum_kg(sys.vA_oovo,t1a,'mnei,ej->mnij')...
% 	     	+0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efij->mnij')...
% 	     	+einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ei->mnif'),t1a,'mnif,fj->mnij');
% 
%      % h2A(abef)
%      H2A.vvvv = sys.vA_vvvv...
% 	     	-einsum_kg(sys.vA_ovvv,t1a,'mbef,am->abef')...
% 	     			+einsum_kg(sys.vA_ovvv,t1a,'maef,bm->abef')...
% 	     	+0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,abmn->abef')...
% 	     	+einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bn->mbef'),t1a,'mbef,am->abef');
%         
%      % h2A(amie)
%      H2A.voov = sys.vA_voov ...
%                 -einsum_kg(sys.vA_ooov,t1a,'nmie,an->amie')...
%                 +einsum_kg(sys.vA_vovv,t1a,'amfe,fi->amie')...
%                 -einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie')...
%                 +einsum_kg(sys.vA_oovv,t2a,'nmfe,afin->amie')...
%                 +einsum_kg(sys.vB_oovv,t2b,'mnef,afin->amie');
%         
%      % h2B(amij)
%      H2B.vooo = sys.vB_vooo...
%                -einsum_kg(sys.vB_oooo,t1a,'nmij,an->amij')...
%                +einsum_kg(sys.vB_vovo,t1a,'amej,ei->amij')...
%                +einsum_kg(sys.vB_voov,t1b,'amie,ej->amij')...
%                +einsum_kg(sys.vB_oovo,t2a,'nmej,aein->amij')...
%                +einsum_kg(sys.vC_oovo,t2b,'nmej,aein->amij')...
%                +einsum_kg(sys.vB_vovv,t2b,'amef,efij->amij')...
%                -einsum_kg(sys.vB_ooov,t2b,'nmif,afnj->amij')...
%                +einsum_kg(sys.fb_ov,t2b,'me,aeij->amij')...
%                -einsum_kg(einsum_kg(sys.vB_oovo,t1a,'nmej,an->amej'),t1a,'amej,ei->amij')...
%                -einsum_kg(einsum_kg(sys.vB_ooov,t1b,'nmie,ej->nmij'),t1a,'nmij,an->amij')...
%                +einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ej->mnjf'),t2b,'mnjf,afin->amij')...
%                +einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t2b,'me,aeij->amij')...
%                -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,an->amfe'),t2b,'amfe,feij->amij')...
%                -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t2b,'nmie,aenj->amij')...
%                +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t2b,'me,aeij->amij')...
%                +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,ej->nmfj'),t2a,'nmfj,afin->amij')...
%                -einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,ej->nmfj'),t1a,'nmfj,an->amfj'),t1a,'amfj,fi->amij')...
%                +einsum_kg(einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie'),t1b,'amie,ej->amij');
% 
%     % h2B(maji)
%     H2B.ovoo =  sys.vB_ovoo...
%                 +einsum_kg(sys.vB_ovvo,t1a,'maei,ej->maji')...
%                 -einsum_kg(sys.vB_oooo,t1b,'mnji,an->maji')...
%                 +einsum_kg(sys.vB_ovov,t1b,'majf,fi->maji')...
%                 +einsum_kg(sys.vB_ooov,t2c,'mnjf,afin->maji')...
%                 +einsum_kg(sys.vB_ovvv,t2b,'maef,efji->maji')...
%                 +einsum_kg(sys.vA_ooov,t2b,'mnjf,fani->maji')...
%                 -einsum_kg(sys.vB_oovo,t2b,'mnei,eajn->maji')...
%                 +einsum_kg(sys.fa_ov,t2b,'me,eaji->maji')...
%                 -einsum_kg(einsum_kg(sys.vB_ooov,t1b,'mnjf,an->majf'),t1b,'majf,fi->maji')...
%                 +einsum_kg(einsum_kg(sys.vB_ovvv,t1a,'maef,ej->majf'),t1b,'majf,fi->maji')...
%                 -einsum_kg(einsum_kg(sys.vB_oovo,t1a,'mnei,ej->mnji'),t1b,'mnji,an->maji')...
%                 +einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ej->mnjf'),t2c,'mnjf,afin->maji')...
%                 +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t2b,'me,eaji->maji')...
%                 +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ej->mnjf'),t2b,'mnjf,fani->maji')...
%                 -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,an->maef'),t2b,'maef,efji->maji')...
%                 +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t2b,'me,eaji->maji')...
%                 -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t2b,'mnei,eajn->maji')...
%                 -einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ej->mnjf'),t1b,'mnjf,an->majf'),t1b,'majf,fi->maji');
% 
%     % h2B(abie)
%     H2B.vvov =  sys.vB_vvov...
% 	    	-einsum_kg(sys.vB_ovov,t1a,'nbie,an->abie')...
% 	    	+einsum_kg(sys.vB_vvvv,t1a,'abfe,fi->abie')...
% 	    	-einsum_kg(sys.vB_voov,t1b,'amie,bm->abie')...
% 	    	+einsum_kg(sys.vB_ovvv,t2a,'nbfe,afin->abie')...
% 	    	+einsum_kg(sys.vC_ovvv,t2b,'nbfe,afin->abie')...
% 	    	+einsum_kg(sys.vB_ooov,t2b,'nmie,abnm->abie')...
% 	    	-einsum_kg(sys.vB_vovv,t2b,'amfe,fbim->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vB_ovvv,t1a,'nbfe,an->abfe'),t1a,'abfe,fi->abie')...
% 	    	+einsum_kg(einsum_kg(sys.vB_ooov,t1a,'nmie,an->amie'),t1b,'amie,bm->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie'),t1b,'amie,bm->abie')...
% 	    	-einsum_kg(sys.fb_ov,t2b,'me,abim->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bm->bnef'),t2b,'bnef,afin->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me'),t2b,'me,abim->abie')...
% 	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t2b,'nmie,abnm->abie')...
% 	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,an->amfe'),t2b,'amfe,fbim->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me'),t2b,'me,abim->abie')...
% 	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmfe,bm->nbfe'),t2a,'nbfe,afin->abie')...
% 	    	+einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,an->amfe'),t1a,'amfe,fi->amie'),t1b,'amie,bm->abie');
% 
%     % h2B(baei)
%     H2B.vvvo =  sys.vB_vvvo...
% 	    	+einsum_kg(sys.vB_vvvv,t1b,'baef,fi->baei')...
% 	    	-einsum_kg(sys.vB_vovo,t1b,'bnei,an->baei')...
% 	    	-einsum_kg(sys.vB_ovvo,t1a,'maei,bm->baei')...
% 	    	+einsum_kg(sys.vB_vovv,t2c,'bnef,afin->baei')...
% 	    	+einsum_kg(sys.vB_oovo,t2b,'mnei,bamn->baei')...
% 	    	+einsum_kg(sys.vA_vovv,t2b,'bnef,fani->baei')...
% 	    	-einsum_kg(sys.vB_ovvv,t2b,'maef,bfmi->baei')...
% 	    	-einsum_kg(sys.fa_ov,t2b,'me,bami->baei')...
% 	    	-einsum_kg(einsum_kg(sys.vB_vovv,t1b,'bnef,fi->bnei'),t1b,'bnei,an->baei')...
% 	    	+einsum_kg(einsum_kg(sys.vB_oovo,t1b,'mnei,an->maei'),t1a,'maei,bm->baei')...
% 	    	-einsum_kg(einsum_kg(sys.vB_ovvv,t1b,'maef,fi->maei'),t1a,'maei,bm->baei')...
% 	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,bm->bnef'),t2c,'bnef,afin->baei')...
% 	    	-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t2b,'me,bami->baei')...
% 	    	-einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bm->bnef'),t2b,'bnef,fani->baei')...
% 	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t2b,'mnei,bamn->baei')...
% 	    	-einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t2b,'me,bami->baei')...
% 	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,an->maef'),t2b,'maef,bfmi->baei')...
% 	    	+einsum_kg(einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,an->maef'),t1b,'maef,fi->maei'),t1a,'maei,bm->baei');
% 
%     % h2B(mnij)
%     H2B.oooo =  sys.vB_oooo...
% 	    	+einsum_kg(sys.vB_oovo,t1a,'mnej,ei->mnij')...
% 	   	+einsum_kg(sys.vB_ooov,t1b,'mnif,fj->mnij')...
% 	   	+einsum_kg(sys.vB_oovv,t2b,'mnef,efij->mnij')...
% 	   	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,ei->mnif'),t1b,'mnif,fj->mnij'); 
% 
%     % h2B(abef)
%     H2B.vvvv =  sys.vB_vvvv...
% 	    	-einsum_kg(sys.vB_ovvv,t1a,'mbef,am->abef')...
% 	    	-einsum_kg(sys.vB_vovv,t1b,'anef,bn->abef')...
% 	    	+einsum_kg(sys.vB_oovv,t2b,'mnef,abmn->abef')...
% 	    	+einsum_kg(einsum_kg(sys.vB_oovv,t1a,'mnef,am->anef'),t1b,'anef,bn->abef');
%         
%     % h2B(amie)
%     H2B.voov = sys.vB_voov ...
%                -einsum_kg(sys.vB_ooov,t1a,'nmie,an->amie')...
%                +einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie')...
%                -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie')...
%                +einsum_kg(sys.vB_oovv,t2a,'nmfe,afin->amie')...
%                +einsum_kg(sys.vC_oovv,t2b,'mnef,afin->amie');
%           
%     % h2B(maie)
%     H2B.ovov = sys.vB_ovov ...
%                +einsum_kg(sys.vB_ovvv,t1a,'mafe,fi->maie')...
%                -einsum_kg(sys.vB_ooov,t1b,'mnie,an->maie')...
%                -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnfe,an->mafe'),t1a,'mafe,fi->maie')...
%                -einsum_kg(sys.vB_oovv,t2b,'mnfe,fain->maie');
%            
%     % h2B(amei)
%     H2B.vovo = sys.vB_vovo ...
%                -einsum_kg(sys.vB_oovo,t1a,'nmei,an->amei')...
%                +einsum_kg(sys.vB_vovv,t1b,'amef,fi->amei')...
%                -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'nmef,fi->nmei'),t1a,'nmei,an->amei')...
%                -einsum_kg(sys.vB_oovv,t2b,'nmef,afni->amei');
%            
%     % h2B(maei)
%     H2B.ovvo = sys.vB_ovvo ...
%                +einsum_kg(sys.vB_ovvv,t1b,'maef,fi->maei')...
%                -einsum_kg(sys.vB_oovo,t1b,'mnei,an->maei')...
%                -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->maei')...
%                +einsum_kg(sys.vB_oovv,t2c,'mnef,afin->maei')...
%                +einsum_kg(sys.vA_oovv,t2b,'mnef,fani->maei');
%    
%      % h2C(amie)
%      H2C.voov = sys.vC_voov ...
%                 -einsum_kg(sys.vC_oovo,t1b,'mnei,an->amie')...
%                 +einsum_kg(sys.vC_vovv,t1b,'amfe,fi->amie')...
%                 -einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->amie')...
%                 +einsum_kg(sys.vC_oovv,t2c,'mnef,afin->amie')...
%                 +einsum_kg(sys.vB_oovv,t2b,'nmfe,fani->amie');

     % (Vt3)_C intermediates
     PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
     PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
     
d1 = -einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3C.PPphHH,'mNef,ABfmNJ->ABeJ')...
+einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3C.PPPhHH,'mNeF,AFBmNJ->ABeJ')...
-einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3C.PPpHhH,'Mnef,ABfMnJ->ABeJ')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3C.PPPHhH,'MneF,AFBMnJ->ABeJ')...
-einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3C.PPpHHH,'MNef,ABfMNJ->ABeJ')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.PPPHHH,'MNeF,AFBMNJ->ABeJ');

d2 = -einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3B.PpPhHH,'Mnef,AfBnMJ->ABeJ')...
-einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3B.PPPhHH,'MneF,AFBnMJ->ABeJ')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3B.PpPHHH,'MNef,AfBMNJ->ABeJ')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PPPHHH,'MNeF,AFBMNJ->ABeJ');

Vt3B.PPvH = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3C.PPpHhH,'mnEf,EBfInJ->mBIJ')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3C.pPPHhH,'mneF,eFBInJ->mBIJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3C.PPPHhH,'mnEF,EFBInJ->mBIJ')...
-einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PPpHHH,'mNEf,EBfINJ->mBIJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPPHHH,'mNeF,eFBINJ->mBIJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPPHHH,'mNEF,EFBINJ->mBIJ');

d2 = +einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3B.PpPhHH,'mneF,FeBnIJ->mBIJ')...
-0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3B.PPPhHH,'mnEF,EFBnIJ->mBIJ')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PpPHHH,'mNeF,FeBINJ->mBIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPPHHH,'mNEF,EFBINJ->mBIJ');

Vt3B.oPHH = d1 + d2;
              
d1 = -einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3C.pPphHH,'mNef,eBfiNJ->mBiJ')...
-einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PPphHH,'mNEf,EBfiNJ->mBiJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPPhHH,'mNeF,eFBiNJ->mBiJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPPhHH,'mNEF,EFBiNJ->mBiJ');

d2 = -einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PpPhHH,'mNeF,FeBiNJ->mBiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPPhHH,'mNEF,EFBiNJ->mBiJ');

Vt3B.oPhH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3B.PpPhHH,'NMfe,AfBiNM->ABie')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3B.PPPhHH,'NMFe,AFBiNM->ABie');

d2 = -0.5*einsum_kg(sys.vC_oovv(HB,HB,pB,:),T3C.PPphHH,'NMfe,ABfiNM->ABie')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),T3C.PPPhHH,'NMFe,AFBiNM->ABie');

Vt3B.PPhv = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3B.PpPHHh,'Nmfe,AfBINm->ABIe')...
+einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3B.PPPHHh,'NmFe,AFBINm->ABIe')...
-einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3B.PpPhHH,'nMfe,AfBnIM->ABIe')...
-einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3B.PPPhHH,'nMFe,AFBnIM->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3B.PpPHHH,'NMfe,AfBINM->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3B.PPPHHH,'NMFe,AFBINM->ABIe');

d2 = -einsum_kg(sys.vC_oovv(hB,HB,pB,:),T3C.PPpHhH,'nMfe,ABfInM->ABIe')...
+einsum_kg(sys.vC_oovv(hB,HB,PB,:),T3C.PPPHhH,'nMFe,AFBInM->ABIe')...
-0.5*einsum_kg(sys.vC_oovv(HB,HB,pB,:),T3C.PPpHHH,'NMfe,ABfINM->ABIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),T3C.PPPHHH,'NMFe,AFBINM->ABIe');

Vt3B.PPHv = -d1 - d2;

d1 =+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3B.PpPhHH,'NmfE,AfEiNJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3B.PPphHH,'NmFe,AFeiNJ->AmiJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3B.PPPhHH,'NmFE,AFEiNJ->AmiJ');

d2 = +einsum_kg(sys.vC_oovv(HB,:,PB,pB),T3C.PPphHH,'NmFe,AFeiNJ->AmiJ')...
+0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),T3C.PPPhHH,'NmFE,AFEiNJ->AmiJ');

Vt3B.PohH = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3B.PpPhHH,'nmfE,AfEnIJ->AmIJ')...
-einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3B.PPphHH,'nmFe,AFenIJ->AmIJ')...
-einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3B.PPPhHH,'nmFE,AFEnIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3B.PpPHHH,'NmfE,AfEINJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3B.PPpHHH,'NmFe,AFeINJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3B.PPPHHH,'NmFE,AFEINJ->AmIJ');

d2 = +einsum_kg(sys.vC_oovv(hB,:,PB,pB),T3C.PPpHhH,'nmFe,AFeInJ->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(hB,:,PB,PB),T3C.PPPHhH,'nmFE,AFEInJ->AmIJ')...
+einsum_kg(sys.vC_oovv(HB,:,PB,pB),T3C.PPpHHH,'NmFe,AFeINJ->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),T3C.PPPHHH,'NmFE,AFEINJ->AmIJ');

Vt3B.PoHH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3B.PPphHH,'MNef,ABfiMN->ABie')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3B.PPPhHH,'MNeF,ABFiMN->ABie');

d2 = +0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3A.PPphHH,'MNef,ABfiMN->ABie')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3A.PPPhHH,'MNeF,ABFiMN->ABie');

Vt3A.PPhv = -d1 - d2;

d1 =-einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3B.PPphHH,'mNef,ABfmIN->ABIe')...
-einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3B.PPPhHH,'mNeF,ABFmIN->ABIe')...
+einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3B.PPpHHh,'Mnef,ABfIMn->ABIe')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3B.PPPHHh,'MneF,ABFIMn->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3B.PPpHHH,'MNef,ABfIMN->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3B.PPPHHH,'MNeF,ABFIMN->ABIe');

d2 = +einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3A.PPphHH,'Mnef,ABfnIM->ABIe')...
+einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3A.PPPhHH,'MneF,ABFnIM->ABIe')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3A.PPpHHH,'MNef,ABfIMN->ABIe')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3A.PPPHHH,'MNeF,ABFIMN->ABIe');

Vt3A.PPHv = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3B.PPphHH,'mNEf,AEfiJN->AmiJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3B.PpPhHH,'mNeF,AeFiJN->AmiJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3B.PPPhHH,'mNEF,AEFiJN->AmiJ');

d2 = -einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3A.PPphHH,'mNeF,AFeiJN->AmiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3A.PPPhHH,'mNEF,AEFiJN->AmiJ');

Vt3A.PohH = d1 + d2;
end