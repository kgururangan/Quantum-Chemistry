function [t3b] = update_t3b_proj6_ccsdt2(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,HBar_t,sys,shift)

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
    
    M23_D1 = einsum_kg(H2B.vvvo(PA,pB,:,HB) + Vt3B.PpvH,t2a(PA,:,hA,HA),'BceK,AeiJ->ABciJK'); % (ab)
    
    M23_D2 = -einsum_kg(I2B_ovoo(:,pB,HA,HB) + Vt3B.opHH,t2a(PA,PA,hA,:),'mcJK,ABim->ABciJK');
    M23_D3 = +einsum_kg(I2B_ovoo(:,pB,hA,HB) + Vt3B.ophH,t2a(PA,PA,HA,:),'mciK,ABJm->ABciJK');
   
    M23_D4 = +einsum_kg(H2B.vvov(PA,pB,hA,:) + Vt3B.Pphv,t2b(PA,:,HA,HB),'Acie,BeJK->ABciJK'); % (ab)
    M23_D5 = -einsum_kg(H2B.vvov(PA,pB,HA,:) + Vt3B.PpHv,t2b(PA,:,hA,HB),'AcJe,BeiK->ABciJK'); % (ab)
    
    M23_D6 = -einsum_kg(I2B_vooo(PA,:,hA,HB) + Vt3B.PohH,t2b(PA,pB,HA,:),'AmiK,BcJm->ABciJK'); % (ab)
    M23_D7 = +einsum_kg(I2B_vooo(PA,:,HA,HB) + Vt3B.PoHH,t2b(PA,pB,hA,:),'AmJK,Bcim->ABciJK'); % (ab)
    
    M23_D8 = +einsum_kg(H2A.vvov(PA,PA,hA,:) + Vt3A.PPhv,t2b(:,pB,HA,HB),'ABie,ecJK->ABciJK');
    M23_D9 = -einsum_kg(H2A.vvov(PA,PA,HA,:) + Vt3A.PPHv,t2b(:,pB,hA,HB),'ABJe,eciK->ABciJK');
    
    M23_D10 = -einsum_kg(I2A_vooo(PA,:,hA,HA) + Vt3A.PohH,t2b(PA,pB,:,HB),'AmiJ,BcmK->ABciJK'); % (ab)
    
    M23_D1456710 = M23_D1 + M23_D4 + M23_D5 + M23_D6 + M23_D7 + M23_D10;
    M23_D1456710 = M23_D1456710 - permute(M23_D1456710,[2,1,3,4,5,6]);
    
    MM23B = M23_D1456710 + M23_D2 + M23_D3 + M23_D8 + M23_D9;
 
    %%%%%
D1 = -einsum_kg(H1A.oo(hA,hA),T3B.PPphHH,'mi,ABcmJK->ABciJK')...
-einsum_kg(H1A.oo(HA,hA),T3B.PPpHHH,'Mi,ABcMJK->ABciJK')...
-einsum_kg(H1A.oo(hA,HA),T3B.PPphhH,'mJ,ABcimK->ABciJK')...
-einsum_kg(H1A.oo(HA,HA),T3B.PPphHH,'MJ,ABciMK->ABciJK');

D2 = -einsum_kg(H1B.oo(hB,HB),T3B.PPphHh,'mK,ABciJm->ABciJK')...
-einsum_kg(H1B.oo(HB,HB),T3B.PPphHH,'MK,ABciJM->ABciJK');
     
D3 = -einsum_kg(H1A.vv(PA,pA),T3B.PpphHH,'Ae,BeciJK->ABciJK')...
+einsum_kg(H1A.vv(PA,PA),T3B.PPphHH,'AE,EBciJK->ABciJK');
D3 = D3 - permute(D3,[2,1,3,4,5,6]);

D4 = +einsum_kg(H1B.vv(pB,pB),T3B.PPphHH,'ce,ABeiJK->ABciJK')...
+einsum_kg(H1B.vv(pB,PB),T3B.PPPhHH,'cE,ABEiJK->ABciJK');

D5 = +0.5*einsum_kg(H2A.oooo(hA,hA,hA,HA),T3B.PPphhH,'mniJ,ABcmnK->ABciJK')...
-einsum_kg(H2A.oooo(HA,hA,hA,HA),T3B.PPphHH,'MniJ,ABcnMK->ABciJK')...
+0.5*einsum_kg(H2A.oooo(HA,HA,hA,HA),T3B.PPpHHH,'MNiJ,ABcMNK->ABciJK');

D6 = +einsum_kg(H2B.oooo(hA,hB,HA,HB),T3B.PPphhh,'mnJK,ABcimn->ABciJK')...
+einsum_kg(H2B.oooo(hA,HB,HA,HB),T3B.PPphhH,'mNJK,ABcimN->ABciJK')...
+einsum_kg(H2B.oooo(HA,hB,HA,HB),T3B.PPphHh,'MnJK,ABciMn->ABciJK')...
+einsum_kg(H2B.oooo(HA,HB,HA,HB),T3B.PPphHH,'MNJK,ABciMN->ABciJK')...
+einsum_kg(H2B.oooo(hA,hB,hA,HB),T3B.PPphHh,'mniK,ABcmJn->ABciJK')...
+einsum_kg(H2B.oooo(hA,HB,hA,HB),T3B.PPphHH,'mNiK,ABcmJN->ABciJK')...
-einsum_kg(H2B.oooo(HA,hB,hA,HB),T3B.PPpHHh,'MniK,ABcJMn->ABciJK')...
-einsum_kg(H2B.oooo(HA,HB,hA,HB),T3B.PPpHHH,'MNiK,ABcJMN->ABciJK');

D7 = +0.5*einsum_kg(H2A.vvvv(PA,PA,pA,pA),T3B.ppphHH,'ABef,efciJK->ABciJK')...
+einsum_kg(H2A.vvvv(PA,PA,PA,pA),T3B.PpphHH,'ABEf,EfciJK->ABciJK')...
+0.5*einsum_kg(H2A.vvvv(PA,PA,PA,PA),T3B.PPphHH,'ABEF,EFciJK->ABciJK');

D8 = +einsum_kg(H2B.vvvv(PA,pB,pA,pB),T3B.PpphHH,'Bcef,AefiJK->ABciJK')...
+einsum_kg(H2B.vvvv(PA,pB,pA,PB),T3B.PpPhHH,'BceF,AeFiJK->ABciJK')...
+einsum_kg(H2B.vvvv(PA,pB,PA,pB),T3B.PPphHH,'BcEf,AEfiJK->ABciJK')...
+einsum_kg(H2B.vvvv(PA,pB,PA,PB),T3B.PPPhHH,'BcEF,AEFiJK->ABciJK');
D8 = D8 - permute(D8,[2,1,3,4,5,6]);
     
D9 = -einsum_kg(H2A.voov(PA,hA,hA,pA),T3B.PpphHH,'Amie,BecmJK->ABciJK')...
+einsum_kg(H2A.voov(PA,hA,hA,PA),T3B.PPphHH,'AmiE,EBcmJK->ABciJK')...
-einsum_kg(H2A.voov(PA,HA,hA,pA),T3B.PppHHH,'AMie,BecMJK->ABciJK')...
+einsum_kg(H2A.voov(PA,HA,hA,PA),T3B.PPpHHH,'AMiE,EBcMJK->ABciJK')...
-einsum_kg(H2A.voov(PA,hA,HA,pA),T3B.PpphhH,'AmJe,BecimK->ABciJK')...
+einsum_kg(H2A.voov(PA,hA,HA,PA),T3B.PPphhH,'AmJE,EBcimK->ABciJK')...
-einsum_kg(H2A.voov(PA,HA,HA,pA),T3B.PpphHH,'AMJe,BeciMK->ABciJK')...
+einsum_kg(H2A.voov(PA,HA,HA,PA),T3B.PPphHH,'AMJE,EBciMK->ABciJK');
    D9 = D9 - permute(D9,[2,1,3,4,5,6]);
    
D10 = -einsum_kg(H2B.voov(PA,hB,hA,pB),T3C.PppHhH,'Amie,BceJmK->ABciJK')...
+einsum_kg(H2B.voov(PA,hB,hA,PB),T3C.PPpHhH,'AmiE,BEcJmK->ABciJK')...
+einsum_kg(H2B.voov(PA,HB,hA,pB),T3C.PppHHH,'AMie,BceJKM->ABciJK')...
-einsum_kg(H2B.voov(PA,HB,hA,PB),T3C.PPpHHH,'AMiE,BEcJKM->ABciJK')...
+einsum_kg(H2B.voov(PA,hB,HA,pB),T3C.PpphhH,'AmJe,BceimK->ABciJK')...
-einsum_kg(H2B.voov(PA,hB,HA,PB),T3C.PPphhH,'AmJE,BEcimK->ABciJK')...
-einsum_kg(H2B.voov(PA,HB,HA,pB),T3C.PpphHH,'AMJe,BceiKM->ABciJK')...
+einsum_kg(H2B.voov(PA,HB,HA,PB),T3C.PPphHH,'AMJE,BEciKM->ABciJK');
    D10 = D10 - permute(D10,[2,1,3,4,5,6]);
    
    D11 = -einsum_kg(H2B.ovvo(hA,pB,pA,HB),T3A.PPphhH,'mceK,ABeimJ->ABciJK')...
-einsum_kg(H2B.ovvo(hA,pB,PA,HB),T3A.PPPhhH,'mcEK,ABEimJ->ABciJK')...
+einsum_kg(H2B.ovvo(HA,pB,pA,HB),T3A.PPphHH,'MceK,ABeiJM->ABciJK')...
+einsum_kg(H2B.ovvo(HA,pB,PA,HB),T3A.PPPhHH,'McEK,ABEiJM->ABciJK');

D12 = +einsum_kg(H2C.voov(pB,hB,HB,pB),T3B.PPphHh,'cmKe,ABeiJm->ABciJK')...
+einsum_kg(H2C.voov(pB,hB,HB,PB),T3B.PPPhHh,'cmKE,ABEiJm->ABciJK')...
+einsum_kg(H2C.voov(pB,HB,HB,pB),T3B.PPphHH,'cMKe,ABeiJM->ABciJK')...
+einsum_kg(H2C.voov(pB,HB,HB,PB),T3B.PPPhHH,'cMKE,ABEiJM->ABciJK');

D13 = +einsum_kg(H2B.vovo(PA,hB,pA,HB),T3B.PpphHh,'AmeK,BeciJm->ABciJK')...
-einsum_kg(H2B.vovo(PA,hB,PA,HB),T3B.PPphHh,'AmEK,EBciJm->ABciJK')...
+einsum_kg(H2B.vovo(PA,HB,pA,HB),T3B.PpphHH,'AMeK,BeciJM->ABciJK')...
-einsum_kg(H2B.vovo(PA,HB,PA,HB),T3B.PPphHH,'AMEK,EBciJM->ABciJK');
D13 = D13 - permute(D13,[2,1,3,4,5,6]);

D14 = -einsum_kg(H2B.ovov(hA,pB,hA,pB),T3B.PPphHH,'mcie,ABemJK->ABciJK')...
-einsum_kg(H2B.ovov(hA,pB,hA,PB),T3B.PPPhHH,'mciE,ABEmJK->ABciJK')...
-einsum_kg(H2B.ovov(HA,pB,hA,pB),T3B.PPpHHH,'Mcie,ABeMJK->ABciJK')...
-einsum_kg(H2B.ovov(HA,pB,hA,PB),T3B.PPPHHH,'MciE,ABEMJK->ABciJK')...
+einsum_kg(H2B.ovov(hA,pB,HA,pB),T3B.PPphhH,'mcJe,ABemiK->ABciJK')...
+einsum_kg(H2B.ovov(hA,pB,HA,PB),T3B.PPPhhH,'mcJE,ABEmiK->ABciJK')...
-einsum_kg(H2B.ovov(HA,pB,HA,pB),T3B.PPphHH,'McJe,ABeiMK->ABciJK')...
-einsum_kg(H2B.ovov(HA,pB,HA,PB),T3B.PPPhHH,'McJE,ABEiMK->ABciJK');

    X3B_PPphHH = MM23B + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14;

    t3b = T3B.PPphHH;
    for a = 1:sys.Nact_p_alpha
        for b = a+1:sys.Nact_p_alpha
            for c = 1:sys.Nunact_p_beta
                for i = 1:sys.Nunact_h_alpha
                    for j = 1:sys.Nact_h_alpha
                        for k = 1:sys.Nact_h_beta
                            denom=sys.fa_hh(i,i)+sys.fa_HH(j,j)+sys.fb_HH(k,k)-sys.fa_PP(a,a)-sys.fa_PP(b,b)-sys.fb_pp(c,c);
                            t3b(a,b,c,i,j,k) = t3b(a,b,c,i,j,k) + X3B_PPphHH(a,b,c,i,j,k)/(denom-shift);
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

d1 = +einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3C.PpphhH,'mnef,AfbmnJ->AbeJ')...
+einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3C.PPphhH,'mneF,AFbmnJ->AbeJ')...
+einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3C.PpphHH,'mNef,AfbmNJ->AbeJ')...
+einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3C.PPphHH,'mNeF,AFbmNJ->AbeJ')...
+einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3C.PppHhH,'Mnef,AfbMnJ->AbeJ')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3C.PPpHhH,'MneF,AFbMnJ->AbeJ')...
+einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3C.PppHHH,'MNef,AfbMNJ->AbeJ')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.PPpHHH,'MNeF,AFbMNJ->AbeJ');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3B.PpphhH,'mnef,AfbmnJ->AbeJ')...
+0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3B.PPphhH,'mneF,AFbmnJ->AbeJ')...
-einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3B.PpphHH,'Mnef,AfbnMJ->AbeJ')...
-einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3B.PPphHH,'MneF,AFbnMJ->AbeJ')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3B.PppHHH,'MNef,AfbMNJ->AbeJ')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PPpHHH,'MNeF,AFbMNJ->AbeJ');

Vt3B.PpvH = -d1 - d2;
              
d1 = +einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3C.pppHhH,'mnef,efbInJ->mbIJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3C.PppHhH,'mnEf,EfbInJ->mbIJ')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3C.pPpHhH,'mneF,eFbInJ->mbIJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3C.PPpHhH,'mnEF,EFbInJ->mbIJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3C.pppHHH,'mNef,efbINJ->mbIJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PppHHH,'mNEf,EfbINJ->mbIJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPpHHH,'mNeF,eFbINJ->mbIJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPpHHH,'mNEF,EFbINJ->mbIJ');

d2 = -0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3B.ppphHH,'mnef,efbnIJ->mbIJ')...
+einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3B.PpphHH,'mneF,FebnIJ->mbIJ')...
-0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3B.PPphHH,'mnEF,EFbnIJ->mbIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3B.pppHHH,'mNef,efbINJ->mbIJ')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PppHHH,'mNeF,FebINJ->mbIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPpHHH,'mNEF,EFbINJ->mbIJ');

Vt3B.opHH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3C.ppphhH,'mnef,efbinJ->mbiJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3C.PpphhH,'mnEf,EfbinJ->mbiJ')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3C.pPphhH,'mneF,eFbinJ->mbiJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3C.PPphhH,'mnEF,EFbinJ->mbiJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3C.ppphHH,'mNef,efbiNJ->mbiJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PpphHH,'mNEf,EfbiNJ->mbiJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPphHH,'mNeF,eFbiNJ->mbiJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPphHH,'mNEF,EFbiNJ->mbiJ');

d2 = +0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3B.ppphhH,'mnef,efbinJ->mbiJ')...
-einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3B.PpphhH,'mneF,FebinJ->mbiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3B.PPphhH,'mnEF,EFbinJ->mbiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3B.ppphHH,'mNef,efbiNJ->mbiJ')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PpphHH,'mNeF,FebiNJ->mbiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPphHH,'mNEF,EFbiNJ->mbiJ');

Vt3B.ophH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3B.Ppphhh,'nmfe,Afbinm->Abie')...
+einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3B.PPphhh,'nmFe,AFbinm->Abie')...
+einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3B.PpphHh,'Nmfe,AfbiNm->Abie')...
+einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3B.PPphHh,'NmFe,AFbiNm->Abie')...
+einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3B.PpphhH,'nMfe,AfbinM->Abie')...
+einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3B.PPphhH,'nMFe,AFbinM->Abie')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3B.PpphHH,'NMfe,AfbiNM->Abie')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3B.PPphHH,'NMFe,AFbiNM->Abie');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,hB,pB,:),T3C.Ppphhh,'nmfe,Afbinm->Abie')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,PB,:),T3C.PPphhh,'nmFe,AFbinm->Abie')...
+einsum_kg(sys.vC_oovv(hB,HB,pB,:),T3C.PpphhH,'nMfe,AfbinM->Abie')...
+einsum_kg(sys.vC_oovv(hB,HB,PB,:),T3C.PPphhH,'nMFe,AFbinM->Abie')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,pB,:),T3C.PpphHH,'NMfe,AfbiNM->Abie')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),T3C.PPphHH,'NMFe,AFbiNM->Abie');

Vt3B.Pphv = -d1 - d2;

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

d1 = -einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3B.PpphHH,'nmfe,AfenIJ->AmIJ')...
-einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3B.PpPhHH,'nmfE,AfEnIJ->AmIJ')...
-einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3B.PPphHH,'nmFe,AFenIJ->AmIJ')...
-einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3B.PPPhHH,'nmFE,AFEnIJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3B.PppHHH,'Nmfe,AfeINJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3B.PpPHHH,'NmfE,AfEINJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3B.PPpHHH,'NmFe,AFeINJ->AmIJ')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3B.PPPHHH,'NmFE,AFEINJ->AmIJ');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,:,pB,pB),T3C.PppHhH,'nmfe,AfeInJ->AmIJ')...
+einsum_kg(sys.vC_oovv(hB,:,PB,pB),T3C.PPpHhH,'nmFe,AFeInJ->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(hB,:,PB,PB),T3C.PPPHhH,'nmFE,AFEInJ->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(HB,:,pB,pB),T3C.PppHHH,'Nmfe,AfeINJ->AmIJ')...
+einsum_kg(sys.vC_oovv(HB,:,PB,pB),T3C.PPpHHH,'NmFe,AFeINJ->AmIJ')...
+0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),T3C.PPPHHH,'NmFE,AFEINJ->AmIJ');

Vt3B.PoHH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3B.PPphhh,'mnef,ABfimn->ABie')...
+einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3B.PPPhhh,'mneF,ABFimn->ABie')...
+einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3B.PPphhH,'mNef,ABfimN->ABie')...
+einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3B.PPPhhH,'mNeF,ABFimN->ABie')...
+einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3B.PPphHh,'Mnef,ABfiMn->ABie')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3B.PPPhHh,'MneF,ABFiMn->ABie')...
+einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3B.PPphHH,'MNef,ABfiMN->ABie')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3B.PPPhHH,'MNeF,ABFiMN->ABie');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3A.PPphhh,'mnef,ABfimn->ABie')...
+0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3A.PPPhhh,'mneF,ABFimn->ABie')...
-einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3A.PPphhH,'Mnef,ABfinM->ABie')...
-einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3A.PPPhhH,'MneF,ABFinM->ABie')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3A.PPphHH,'MNef,ABfiMN->ABie')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3A.PPPhHH,'MNeF,ABFiMN->ABie');

Vt3A.PPhv = -d1 - d2;

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

d1 = +einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3B.PpphHh,'mnef,AefiJn->AmiJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3B.PPphHh,'mnEf,AEfiJn->AmiJ')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3B.PpPhHh,'mneF,AeFiJn->AmiJ')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3B.PPPhHh,'mnEF,AEFiJn->AmiJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3B.PpphHH,'mNef,AefiJN->AmiJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3B.PPphHH,'mNEf,AEfiJN->AmiJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3B.PpPhHH,'mNeF,AeFiJN->AmiJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3B.PPPhHH,'mNEF,AEFiJN->AmiJ');

d2 = -0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3A.PpphhH,'mnef,AefinJ->AmiJ')...
+einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3A.PPphhH,'mneF,AFeinJ->AmiJ')...
-0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3A.PPPhhH,'mnEF,AEFinJ->AmiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3A.PpphHH,'mNef,AefiJN->AmiJ')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3A.PPphHH,'mNeF,AFeiJN->AmiJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3A.PPPhHH,'mNEF,AEFiJN->AmiJ');

Vt3A.PohH = d1 + d2;

end