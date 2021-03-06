function [t3c] = update_t3c_proj6_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift)
    
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    [Vt3B,Vt3C] = get_t3c_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys);
    
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
    
    % MM23C
    I2B_vooo = H2B.vooo - einsum_kg(H1B.ov,t2b,'me,aeij->amij');
    I2C_ovoo = H2C.ovoo - einsum_kg(H1B.ov,t2c,'me,ecjk->mcjk');
    I2B_ovoo = H2B.ovoo - einsum_kg(H1A.ov,t2b,'me,ebij->mbij');
    

    M23_D1 = +einsum_kg(H2B.vvov(PA,PB,HA,:) + Vt3B.PPHv,t2c(:,pB,hB,HB),'ABIe,ecjK->ABcIjK');
    M23_D2 = -einsum_kg(H2B.vvov(PA,pB,HA,:) + Vt3B.PpHv,t2c(:,PB,hB,HB),'AcIe,eBjK->ABcIjK');
    
    M23_D3 = -einsum_kg(I2B_vooo(PA,:,HA,hB) + Vt3B.PoHh,t2c(PB,pB,:,HB),'AmIj,BcmK->ABcIjK');
    M23_D4 = +einsum_kg(I2B_vooo(PA,:,HA,HB) + Vt3B.PoHH,t2c(PB,pB,:,hB),'AmIK,Bcmj->ABcIjK');
  
    M23_D5 = +einsum_kg(H2C.vvvo(PB,pB,:,HB) + Vt3C.PpvH,t2b(PA,:,HA,hB),'BceK,AeIj->ABcIjK');
    M23_D6 = -einsum_kg(H2C.vvvo(PB,pB,:,hB) + Vt3C.Ppvh,t2b(PA,:,HA,HB),'Bcej,AeIK->ABcIjK');

    M23_D7 = -einsum_kg(I2C_ovoo(:,pB,hB,HB) + Vt3C.ophH,t2b(PA,PB,HA,:),'mcjK,ABIm->ABcIjK');
    M23_D8 = +einsum_kg(I2C_ovoo(:,PB,hB,HB) + Vt3C.oPhH,t2b(PA,pB,HA,:),'mBjK,AcIm->ABcIjK');

    M23_D9 = +einsum_kg(H2B.vvvo(PA,PB,:,hB) + Vt3B.PPvh,t2b(:,pB,HA,HB),'ABej,ecIK->ABcIjK');
    M23_D10 = -einsum_kg(H2B.vvvo(PA,pB,:,hB) + Vt3B.Ppvh,t2b(:,PB,HA,HB),'Acej,eBIK->ABcIjK');
    M23_D11 = -einsum_kg(H2B.vvvo(PA,PB,:,HB) + Vt3B.PPvH,t2b(:,pB,HA,hB),'ABeK,ecIj->ABcIjK');
    M23_D12 = +einsum_kg(H2B.vvvo(PA,pB,:,HB) + Vt3B.PpvH,t2b(:,PB,HA,hB),'AceK,eBIj->ABcIjK');

    M23_D13 = -einsum_kg(I2B_ovoo(:,PB,HA,hB) + Vt3B.oPHh,t2b(PA,pB,:,HB),'mBIj,AcmK->ABcIjK');
    M23_D14 = +einsum_kg(I2B_ovoo(:,pB,HA,hB) + Vt3B.opHh,t2b(PA,PB,:,HB),'mcIj,ABmK->ABcIjK');
    M23_D15 = +einsum_kg(I2B_ovoo(:,PB,HA,HB) + Vt3B.oPHH,t2b(PA,pB,:,hB),'mBIK,Acmj->ABcIjK');
    M23_D16 = -einsum_kg(I2B_ovoo(:,pB,HA,HB) + Vt3B.opHH,t2b(PA,PB,:,hB),'mcIK,ABmj->ABcIjK');

    MM23C = M23_D1 + M23_D2 + M23_D3 + M23_D4 + M23_D5 + M23_D6 + M23_D7 + ...
            M23_D8 + M23_D9 + M23_D10 + M23_D11 + M23_D12 + M23_D13 + M23_D14 + ...
            M23_D15 + M23_D16;
        
    % (HBar*T3)_C

D1 = -einsum_kg(H1A.oo(HA,HA),T3C.PPpHhH,'MI,ABcMjK->ABcIjK');

D2 = -einsum_kg(H1B.oo(hB,hB),T3C.PPpHhH,'mj,ABcImK->ABcIjK')...
-einsum_kg(H1B.oo(HB,hB),T3C.PPpHHH,'Mj,ABcIMK->ABcIjK')...
-einsum_kg(H1B.oo(HB,HB),T3C.PPpHhH,'MK,ABcIjM->ABcIjK');

D3 = +einsum_kg(H1A.vv(PA,PA),T3C.PPpHhH,'AE,EBcIjK->ABcIjK');

D4 =+einsum_kg(H1B.vv(PB,PB),T3C.PPpHhH,'BE,AEcIjK->ABcIjK')...
+einsum_kg(H1B.vv(pB,pB),T3C.PPpHhH,'ce,ABeIjK->ABcIjK')...
-einsum_kg(H1B.vv(pB,PB),T3C.PPPHhH,'cE,AEBIjK->ABcIjK');

D5 = -einsum_kg(H2C.oooo(HB,hB,hB,HB),T3C.PPpHhH,'MnjK,ABcInM->ABcIjK')...
+0.5*einsum_kg(H2C.oooo(HB,HB,hB,HB),T3C.PPpHHH,'MNjK,ABcIMN->ABcIjK');

D6 = +einsum_kg(H2B.oooo(hA,HB,HA,hB),T3C.PPphHH,'mNIj,ABcmNK->ABcIjK')...
+einsum_kg(H2B.oooo(HA,hB,HA,hB),T3C.PPpHhH,'MnIj,ABcMnK->ABcIjK')...
+einsum_kg(H2B.oooo(HA,HB,HA,hB),T3C.PPpHHH,'MNIj,ABcMNK->ABcIjK')...
+einsum_kg(H2B.oooo(HA,HB,HA,HB),T3C.PPpHhH,'MNIK,ABcMjN->ABcIjK');

D7 = +einsum_kg(H2C.vvvv(PB,pB,PB,pB),T3C.PPpHhH,'BcEf,AEfIjK->ABcIjK')...
+0.5*einsum_kg(H2C.vvvv(PB,pB,PB,PB),T3C.PPPHhH,'BcEF,AEFIjK->ABcIjK');

D8 = +einsum_kg(H2B.vvvv(PA,PB,PA,PB),T3C.PPpHhH,'ABEF,EFcIjK->ABcIjK')...
-einsum_kg(H2B.vvvv(PA,pB,pA,PB),T3C.pPPHhH,'AceF,eFBIjK->ABcIjK')...
+einsum_kg(H2B.vvvv(PA,pB,PA,pB),T3C.PPpHhH,'AcEf,EBfIjK->ABcIjK')...
-einsum_kg(H2B.vvvv(PA,pB,PA,PB),T3C.PPPHhH,'AcEF,EFBIjK->ABcIjK');

D9 = +einsum_kg(H2A.voov(PA,HA,HA,PA),T3C.PPpHhH,'AMIE,EBcMjK->ABcIjK');

D10 = -einsum_kg(H2B.voov(PA,HB,HA,PB),T3D.PPphHH,'AMIE,EBcjMK->ABcIjK');

D11 = -einsum_kg(H2B.ovvo(hA,PB,PA,hB),T3B.PPphHH,'mBEj,AEcmIK->ABcIjK')...
+einsum_kg(H2B.ovvo(HA,PB,PA,hB),T3B.PPpHHH,'MBEj,AEcIMK->ABcIjK')...
+einsum_kg(H2B.ovvo(hA,pB,pA,hB),T3B.PpPhHH,'mcej,AeBmIK->ABcIjK')...
+einsum_kg(H2B.ovvo(hA,pB,PA,hB),T3B.PPPhHH,'mcEj,AEBmIK->ABcIjK')...
-einsum_kg(H2B.ovvo(HA,pB,pA,hB),T3B.PpPHHH,'Mcej,AeBIMK->ABcIjK')...
-einsum_kg(H2B.ovvo(HA,pB,PA,hB),T3B.PPPHHH,'McEj,AEBIMK->ABcIjK')...
-einsum_kg(H2B.ovvo(HA,PB,PA,HB),T3B.PPpHHh,'MBEK,AEcIMj->ABcIjK')...
+einsum_kg(H2B.ovvo(HA,pB,pA,HB),T3B.PpPHHh,'MceK,AeBIMj->ABcIjK')...
+einsum_kg(H2B.ovvo(HA,pB,PA,HB),T3B.PPPHHh,'McEK,AEBIMj->ABcIjK');

D12 = +einsum_kg(H2C.voov(PB,hB,hB,PB),T3C.PPpHhH,'BmjE,AEcImK->ABcIjK')...
+einsum_kg(H2C.voov(PB,HB,hB,PB),T3C.PPpHHH,'BMjE,AEcIMK->ABcIjK')...
+einsum_kg(H2C.voov(pB,hB,hB,pB),T3C.PPpHhH,'cmje,ABeImK->ABcIjK')...
-einsum_kg(H2C.voov(pB,hB,hB,PB),T3C.PPPHhH,'cmjE,AEBImK->ABcIjK')...
+einsum_kg(H2C.voov(pB,HB,hB,pB),T3C.PPpHHH,'cMje,ABeIMK->ABcIjK')...
-einsum_kg(H2C.voov(pB,HB,hB,PB),T3C.PPPHHH,'cMjE,AEBIMK->ABcIjK')...
+einsum_kg(H2C.voov(PB,HB,HB,PB),T3C.PPpHhH,'BMKE,AEcIjM->ABcIjK')...
+einsum_kg(H2C.voov(pB,HB,HB,pB),T3C.PPpHhH,'cMKe,ABeIjM->ABcIjK')...
-einsum_kg(H2C.voov(pB,HB,HB,PB),T3C.PPPHhH,'cMKE,AEBIjM->ABcIjK');

D13 = -einsum_kg(H2B.ovov(HA,PB,HA,PB),T3C.PPpHhH,'MBIE,AEcMjK->ABcIjK')...
-einsum_kg(H2B.ovov(HA,pB,HA,pB),T3C.PPpHhH,'McIe,ABeMjK->ABcIjK')...
+einsum_kg(H2B.ovov(HA,pB,HA,PB),T3C.PPPHhH,'McIE,AEBMjK->ABcIjK');

D14 = -einsum_kg(H2B.vovo(PA,hB,PA,hB),T3C.PPpHhH,'AmEj,EBcImK->ABcIjK')...
-einsum_kg(H2B.vovo(PA,HB,PA,hB),T3C.PPpHHH,'AMEj,EBcIMK->ABcIjK')...
-einsum_kg(H2B.vovo(PA,HB,PA,HB),T3C.PPpHhH,'AMEK,EBcIjM->ABcIjK');

X3C_PPpHhH =  D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14 + MM23C;

    t3c = T3C.PPpHhH;
    for a = 1:sys.Nact_p_alpha
        for b = 1:sys.Nact_p_beta
            for c = 1:sys.Nunact_p_beta
                for i = 1:sys.Nact_h_alpha
                    for j = 1:sys.Nunact_h_beta
                        for k = 1:sys.Nact_h_beta
                            denom=sys.fa_HH(i,i)+sys.fb_hh(j,j)+sys.fb_HH(k,k)-sys.fa_PP(a,a)-sys.fb_PP(b,b)-sys.fb_pp(c,c);
                            t3c(a,b,c,i,j,k) = t3c(a,b,c,i,j,k) + X3C_PPpHhH(a,b,c,i,j,k)/(denom-shift);
                        end
                    end
                end
            end
        end
    end

end

function [Vt3B,Vt3C] = get_t3c_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys)

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
%     % h2A(amie)
%     H2A.voov = sys.vA_voov ...
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

     % (Vt3)_C intermediates
     PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
     PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
     
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
     
d1 = -einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3B.PPphHh,'nmFe,AFbnIm->AbIe')...
+einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3B.PPpHHh,'NmFe,AFbINm->AbIe')...
-einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3B.PPphHH,'nMFe,AFbnIM->AbIe')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3B.PPpHHH,'NMFe,AFbINM->AbIe');

d2 = +einsum_kg(sys.vC_oovv(hB,HB,PB,:),T3C.PPpHhH,'nMFe,AFbInM->AbIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),T3C.PPpHHH,'NMFe,AFbINM->AbIe');

Vt3B.PpHv = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3B.PpPHHh,'NmfE,AfEINj->AmIj')...
+einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3B.PPpHHh,'NmFe,AFeINj->AmIj')...
+einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3B.PPPHHh,'NmFE,AFEINj->AmIj');

d2 = -einsum_kg(sys.vC_oovv(HB,:,PB,pB),T3C.PPpHhH,'NmFe,AFeIjN->AmIj')...
-0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),T3C.PPPHhH,'NmFE,AFEIjN->AmIj');

Vt3B.PoHh = d1 + d2;

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

d1 = +einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3C.PPpHhH,'NmFe,FBaNmI->BaeI')...
-einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3C.PPphHH,'nMFe,FBanIM->BaeI')...
-einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPpHHH,'NMFe,FBaNIM->BaeI');

d2 = +einsum_kg(sys.vC_oovv(HB,hB,:,PB),T3D.PPphHH,'MneF,BFanIM->BaeI')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPpHHH,'MNeF,BFaIMN->BaeI');

Vt3C.PpvH = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3C.PPpHhH,'NMFe,FBaNiM->Baei');

d2 = +0.5*einsum_kg(sys.vC_oovv(HB,HB,:,PB),T3D.PPphHH,'MNeF,BFaiMN->Baei');

Vt3C.Ppvh = -d1 - d2;


d1 = +einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPpHhH,'NmFE,FEaNjI->majI');

d2 = -0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPphHH,'mNEF,EFajIN->majI');

Vt3C.ophH = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3C.pPPHhH,'NmfE,fAENjI->mAjI')...
-einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3C.PPpHhH,'NmFe,FAeNjI->mAjI')...
-einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3C.PPPHhH,'NmFE,FAENjI->mAjI');

d2 = +einsum_kg(sys.vC_oovv(:,HB,pB,PB),T3D.PPphHH,'mNeF,AFejIN->mAjI')...
-0.5*einsum_kg(sys.vC_oovv(:,HB,PB,PB),T3D.PPPhHH,'mNEF,AEFjIN->mAjI');

Vt3C.oPhH = d1 + d2;

d1 = +einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3C.PPpHhH,'MNef,ABfMjN->ABej')...
-einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.PPPHhH,'MNeF,AFBMjN->ABej');

d2 = +0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3B.PpPHHh,'MNef,AfBMNj->ABej')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PPPHHh,'MNeF,AFBMNj->ABej');

Vt3B.PPvh = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.PPpHhH,'MNeF,AFbMjN->Abej');

d2 = +0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PPpHHh,'MNeF,AFbMNj->Abej');

Vt3B.Ppvh = -d1 - d2;

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

d1 = +einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3C.PPphHH,'mNeF,AFbmNJ->AbeJ')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3C.PPpHhH,'MneF,AFbMnJ->AbeJ')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.PPpHHH,'MNeF,AFbMNJ->AbeJ');

d2 = -einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3B.PPphHH,'MneF,AFbnMJ->AbeJ')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PPpHHH,'MNeF,AFbMNJ->AbeJ');

Vt3B.PpvH = -d1 - d2;

d1 = +einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PPpHhH,'mNEf,EBfIjN->mBIj')...
-einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPPHhH,'mNeF,eFBIjN->mBIj')...
-einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPPHhH,'mNEF,EFBIjN->mBIj');

d2 = -einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PpPHHh,'mNeF,FeBINj->mBIj')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPPHHh,'mNEF,EFBINj->mBIj');

Vt3B.oPHh = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPpHhH,'mNEF,EFbIjN->mbIj');

d2 = +0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPpHHh,'mNEF,EFbINj->mbIj');

Vt3B.opHh = d1 + d2;

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

d1 = +einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3C.PPpHhH,'mnEF,EFbInJ->mbIJ')...
+einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPpHHH,'mNEF,EFbINJ->mbIJ');

d2 = -0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3B.PPphHH,'mnEF,EFbnIJ->mbIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPpHHH,'mNEF,EFbINJ->mbIJ');

Vt3B.opHH = d1 + d2;
end