function [t3b] = update_t3b_proj9_ccsdt2(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,HBar_t,sys,shift)

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
    
    M23_D1 = +einsum_kg(H2B.vvvo(pA,PB,:,hB) + Vt3B.pPvh,t2a(PA,:,HA,HA),'bCek,AeIJ->AbCIJk');
    M23_D2 = -einsum_kg(H2B.vvvo(PA,PB,:,hB) + Vt3B.PPvh,t2a(pA,:,HA,HA),'ACek,beIJ->AbCIJk');
    
    M23_D3 = -einsum_kg(I2B_ovoo(:,PB,HA,hB) + Vt3B.oPHh,t2a(PA,pA,HA,:),'mCJk,AbIm->AbCIJk'); % A(ij)
    
    M23_D4 = +einsum_kg(H2B.vvov(PA,PB,HA,:) + Vt3B.PPHv,t2b(pA,:,HA,hB),'ACIe,beJk->AbCIJk'); % A(ij)
    M23_D5 = -einsum_kg(H2B.vvov(pA,PB,HA,:) + Vt3B.pPHv,t2b(PA,:,HA,hB),'bCIe,AeJk->AbCIJk'); % A(ij)
    
    M23_D6 = -einsum_kg(I2B_vooo(PA,:,HA,hB) + Vt3B.PoHh,t2b(pA,PB,HA,:),'AmIk,bCJm->AbCIJk'); % A(ij)
    M23_D7 = +einsum_kg(I2B_vooo(pA,:,HA,hB) + Vt3B.poHh,t2b(PA,PB,HA,:),'bmIk,ACJm->AbCIJk'); % A(ij)
    
    M23_D8 = +einsum_kg(H2A.vvov(PA,pA,HA,:) + Vt3A.PpHv,t2b(:,PB,HA,hB),'AbIe,eCJk->AbCIJk'); % A(ij)
    
    M23_D9 = -einsum_kg(I2A_vooo(PA,:,HA,HA) + Vt3A.PoHH,t2b(pA,PB,:,hB),'AmIJ,bCmk->AbCIJk');
    M23_D10 = +einsum_kg(I2A_vooo(pA,:,HA,HA) + Vt3A.poHH,t2b(PA,PB,:,hB),'bmIJ,ACmk->AbCIJk');
    

    M23_D345678 = M23_D3 + M23_D4 + M23_D5 + M23_D6 + M23_D7 + M23_D8;
    M23_D345678 = M23_D345678 - permute(M23_D345678,[1,2,3,5,4,6]);

    MM23B = M23_D345678 + M23_D1 + M23_D2 + M23_D9 + M23_D10;
    
 
    %%%%%
    
    D1 = -einsum_kg(H1A.oo(hA,HA),T3B.PpPhHh,'mI,AbCmJk->AbCIJk')...
    -einsum_kg(H1A.oo(HA,HA),T3B.PpPHHh,'MI,AbCMJk->AbCIJk');
    D1 = D1 - permute(D1,[1,2,3,5,4,6]);

    D2 = -einsum_kg(H1B.oo(hB,hB),T3B.PpPHHh,'mk,AbCIJm->AbCIJk')...
    -einsum_kg(H1B.oo(HB,hB),T3B.PpPHHH,'Mk,AbCIJM->AbCIJk');

    D3 = +einsum_kg(H1A.vv(PA,pA),T3B.ppPHHh,'Ae,ebCIJk->AbCIJk')...
    +einsum_kg(H1A.vv(PA,PA),T3B.PpPHHh,'AE,EbCIJk->AbCIJk');

    D4 = +einsum_kg(H1A.vv(pA,PA),T3B.PPPHHh,'bE,AECIJk->AbCIJk')...
    +einsum_kg(H1A.vv(pA,pA),T3B.PpPHHh,'be,AeCIJk->AbCIJk');

    D5 = +einsum_kg(H1B.vv(PB,PB),T3B.PpPHHh,'CE,AbEIJk->AbCIJk')...
    +einsum_kg(H1B.vv(PB,pB),T3B.PppHHh,'Ce,AbeIJk->AbCIJk');

    D6 = +einsum_kg(H2A.voov(PA,HA,HA,PA),T3B.PpPHHh,'AMIE,EbCMJk->AbCIJk')...
    +einsum_kg(H2A.voov(PA,hA,HA,PA),T3B.PpPhHh,'AmIE,EbCmJk->AbCIJk')...
    +einsum_kg(H2A.voov(PA,HA,HA,pA),T3B.ppPHHh,'AMIe,ebCMJk->AbCIJk')...
    +einsum_kg(H2A.voov(PA,hA,HA,pA),T3B.ppPhHh,'AmIe,ebCmJk->AbCIJk');

    D7 = -einsum_kg(H2A.voov(pA,hA,HA,pA),T3B.PpPhHh,'bmJe,AeCmIk->AbCIJk')...
    -einsum_kg(H2A.voov(pA,hA,HA,PA),T3B.PPPhHh,'bmJE,AECmIk->AbCIJk')...
    +einsum_kg(H2A.voov(pA,HA,HA,pA),T3B.PpPHHh,'bMJe,AeCIMk->AbCIJk')...
    +einsum_kg(H2A.voov(pA,HA,HA,PA),T3B.PPPHHh,'bMJE,AECIMk->AbCIJk');

    D8 = +einsum_kg(H2B.ovvo(HA,PB,PA,hB),T3A.PPpHHH,'MCEk,AEbIMJ->AbCIJk')...
    -einsum_kg(H2B.ovvo(hA,PB,PA,hB),T3A.PPphHH,'mCEk,AEbmIJ->AbCIJk')...
    +einsum_kg(H2B.ovvo(HA,PB,pA,hB),T3A.PppHHH,'MCek,AbeIJM->AbCIJk')...
    -einsum_kg(H2B.ovvo(hA,PB,pA,hB),T3A.PpphHH,'mCek,AbemJI->AbCIJk');

    D9 = +einsum_kg(H2B.voov(PA,HB,HA,PB),T3C.pPPHhH,'AMIE,bCEJkM->AbCIJk')...
    +einsum_kg(H2B.voov(PA,hB,HA,PB),T3C.pPPHhh,'AmIE,bCEJkm->AbCIJk')...
    +einsum_kg(H2B.voov(PA,HB,HA,pB),T3C.pPpHhH,'AMIe,bCeJkM->AbCIJk')...
    +einsum_kg(H2B.voov(PA,hB,HA,pB),T3C.pPpHhh,'AmIe,bCeJkm->AbCIJk');

    D10 = +einsum_kg(H2B.voov(pA,HB,HA,PB),T3C.PPPHhH,'bMJE,ACEIkM->AbCIJk')...
    +einsum_kg(H2B.voov(pA,hB,HA,PB),T3C.PPPHhh,'bmJE,ACEIkm->AbCIJk')...
    +einsum_kg(H2B.voov(pA,HB,HA,pB),T3C.PPpHhH,'bMJe,ACeIkM->AbCIJk')...
    +einsum_kg(H2B.voov(pA,hB,HA,pB),T3C.PPpHhh,'bmJe,ACeIkm->AbCIJk');

    D11 = +einsum_kg(H2C.voov(PB,HB,hB,PB),T3B.PpPHHH,'CMkE,AbEIJM->AbCIJk')...
    +einsum_kg(H2C.voov(PB,hB,hB,PB),T3B.PpPHHh,'CmkE,AbEIJm->AbCIJk')...
    +einsum_kg(H2C.voov(PB,HB,hB,pB),T3B.PppHHH,'CMke,AbeIJM->AbCIJk')...
    +einsum_kg(H2C.voov(PB,hB,hB,pB),T3B.PppHHh,'Cmke,AbeIJm->AbCIJk');

    D12 = -einsum_kg(H2B.vovo(PA,HB,PA,hB),T3B.PpPHHH,'AMEk,EbCIJM->AbCIJk')...
    -einsum_kg(H2B.vovo(PA,hB,PA,hB),T3B.PpPHHh,'AmEk,EbCIJm->AbCIJk')...
    -einsum_kg(H2B.vovo(PA,HB,pA,hB),T3B.ppPHHH,'AMek,ebCIJM->AbCIJk')...
    -einsum_kg(H2B.vovo(PA,hB,pA,hB),T3B.ppPHHh,'Amek,ebCIJm->AbCIJk');

    D13 = -einsum_kg(H2B.vovo(pA,HB,PA,hB),T3B.PPPHHH,'bMEk,AECIJM->AbCIJk')...
    -einsum_kg(H2B.vovo(pA,hB,PA,hB),T3B.PPPHHh,'bmEk,AECIJm->AbCIJk')...
    -einsum_kg(H2B.vovo(pA,HB,pA,hB),T3B.PpPHHH,'bMek,AeCIJM->AbCIJk')...
    -einsum_kg(H2B.vovo(pA,hB,pA,hB),T3B.PpPHHh,'bmek,AeCIJm->AbCIJk');

    D14 = -einsum_kg(H2B.ovov(HA,PB,HA,PB),T3B.PpPHHh,'MCIE,AbEMJk->AbCIJk')...
    -einsum_kg(H2B.ovov(hA,PB,HA,PB),T3B.PpPhHh,'mCIE,AbEmJk->AbCIJk')...
    -einsum_kg(H2B.ovov(HA,PB,HA,pB),T3B.PppHHh,'MCIe,AbeMJk->AbCIJk')...
    -einsum_kg(H2B.ovov(hA,PB,HA,pB),T3B.PpphHh,'mCIe,AbemJk->AbCIJk');

    D15 = +0.5*einsum_kg(H2A.oooo(HA,HA,HA,HA),T3B.PpPHHh,'MNIJ,AbCMNk->AbCIJk')...
    +einsum_kg(H2A.oooo(hA,HA,HA,HA),T3B.PpPhHh,'mNIJ,AbCmNk->AbCIJk')...
    +0.5*einsum_kg(H2A.oooo(hA,hA,HA,HA),T3B.PpPhhh,'mnIJ,AbCmnk->AbCIJk');

    D16 = +einsum_kg(H2B.oooo(HA,HB,HA,hB),T3B.PpPHHH,'MNIk,AbCMJN->AbCIJk')...
    +einsum_kg(H2B.oooo(hA,HB,HA,hB),T3B.PpPhHH,'mNIk,AbCmJN->AbCIJk')...
    +einsum_kg(H2B.oooo(HA,hB,HA,hB),T3B.PpPHHh,'MnIk,AbCMJn->AbCIJk')...
    +einsum_kg(H2B.oooo(hA,hB,HA,hB),T3B.PpPhHh,'mnIk,AbCmJn->AbCIJk');

    D17 = +0.5*einsum_kg(H2A.vvvv(PA,pA,PA,PA),T3B.PPPHHh,'AbEF,EFCIJk->AbCIJk')...
    +einsum_kg(H2A.vvvv(PA,pA,PA,pA),T3B.PpPHHh,'AbEf,EfCIJk->AbCIJk')...
    +0.5*einsum_kg(H2A.vvvv(PA,pA,pA,pA),T3B.ppPHHh,'Abef,efCIJk->AbCIJk');

    D18 = +einsum_kg(H2B.vvvv(PA,PB,PA,PB),T3B.PpPHHh,'ACEF,EbFIJk->AbCIJk')...
    +einsum_kg(H2B.vvvv(PA,PB,PA,pB),T3B.PppHHh,'ACEf,EbfIJk->AbCIJk')...
    +einsum_kg(H2B.vvvv(PA,PB,pA,PB),T3B.ppPHHh,'ACeF,ebFIJk->AbCIJk')...
    +einsum_kg(H2B.vvvv(PA,PB,pA,pB),T3B.pppHHh,'ACef,ebfIJk->AbCIJk');

    D19 = +einsum_kg(H2B.vvvv(pA,PB,PA,PB),T3B.PPPHHh,'bCEF,AEFIJk->AbCIJk')...
    +einsum_kg(H2B.vvvv(pA,PB,PA,pB),T3B.PPpHHh,'bCEf,AEfIJk->AbCIJk')...
    +einsum_kg(H2B.vvvv(pA,PB,pA,PB),T3B.PpPHHh,'bCeF,AeFIJk->AbCIJk')...
    +einsum_kg(H2B.vvvv(pA,PB,pA,pB),T3B.PppHHh,'bCef,AefIJk->AbCIJk');

    D_6_7_9_10_14_16 = D6 + D7 + D9 + D10 + D14 + D16;
    D_6_7_9_10_14_16 = D_6_7_9_10_14_16 - permute(D_6_7_9_10_14_16,[1,2,3,5,4,6]);

    X3B_PpPHHh = MM23B + D1 + D2 + D3 + D4 + D5 + D_6_7_9_10_14_16 + D8 + D11 + D12 + D13 + D15 + D17 + D18 + D19;

    t3b = T3B.PpPHHh;
    for a = 1:sys.Nact_p_alpha
        for b = 1:sys.Nunact_p_alpha
            for c = 1:sys.Nact_p_beta
                for i = 1:sys.Nact_h_alpha
                    for j = i+1:sys.Nact_h_alpha
                        for k = 1:sys.Nunact_h_beta
                            denom = sys.fa_HH(i,i)+sys.fa_HH(j,j)+sys.fb_hh(k,k)-sys.fa_PP(a,a)-sys.fa_pp(b,b)-sys.fb_PP(c,c);
                            t3b(a,b,c,i,j,k) = t3b(a,b,c,i,j,k) + X3B_PpPHHh(a,b,c,i,j,k)/(denom-shift);
                            t3b(a,b,c,j,i,k) = -t3b(a,b,c,i,j,k);
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
    
d1 = -einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3C.pPphhh,'mnef,aBfmnj->aBej')...
+einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3C.pPPhhh,'mneF,aFBmnj->aBej')...
+einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3C.pPphhH,'mNef,aBfmjN->aBej')...
-einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3C.pPPhhH,'mNeF,aFBmjN->aBej')...
-einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3C.pPpHhh,'Mnef,aBfMnj->aBej')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3C.pPPHhh,'MneF,aFBMnj->aBej')...
+einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3C.pPpHhH,'MNef,aBfMjN->aBej')...
-einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.pPPHhH,'MNeF,aFBMjN->aBej');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3B.ppPhhh,'mnef,afBmnj->aBej')...
-0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3B.PpPhhh,'mneF,FaBmnj->aBej')...
-einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3B.ppPhHh,'Mnef,afBnMj->aBej')...
+einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3B.PpPhHh,'MneF,FaBnMj->aBej')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3B.ppPHHh,'MNef,afBMNj->aBej')...
-0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PpPHHh,'MNeF,FaBMNj->aBej');

Vt3B.pPvh = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3C.PPphhh,'mnef,ABfmnj->ABej')...
+einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3C.PPPhhh,'mneF,AFBmnj->ABej')...
+einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3C.PPphhH,'mNef,ABfmjN->ABej')...
-einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3C.PPPhhH,'mNeF,AFBmjN->ABej')...
-einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3C.PPpHhh,'Mnef,ABfMnj->ABej')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3C.PPPHhh,'MneF,AFBMnj->ABej')...
+einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3C.PPpHhH,'MNef,ABfMjN->ABej')...
-einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3C.PPPHhH,'MNeF,AFBMjN->ABej');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3B.PpPhhh,'mnef,AfBmnj->ABej')...
+0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3B.PPPhhh,'mneF,AFBmnj->ABej')...
-einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3B.PpPhHh,'Mnef,AfBnMj->ABej')...
-einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3B.PPPhHh,'MneF,AFBnMj->ABej')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3B.PpPHHh,'MNef,AfBMNj->ABej')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3B.PPPHHh,'MNeF,AFBMNj->ABej');

Vt3B.PPvh = -d1 - d2;

d1 = -einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3C.pPpHhh,'mnef,eBfInj->mBIj')...
-einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3C.PPpHhh,'mnEf,EBfInj->mBIj')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3C.pPPHhh,'mneF,eFBInj->mBIj')...
+einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3C.PPPHhh,'mnEF,EFBInj->mBIj')...
+einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3C.pPpHhH,'mNef,eBfIjN->mBIj')...
+einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3C.PPpHhH,'mNEf,EBfIjN->mBIj')...
-einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3C.pPPHhH,'mNeF,eFBIjN->mBIj')...
-einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3C.PPPHhH,'mNEF,EFBIjN->mBIj');

d2 = -0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3B.ppPhHh,'mnef,efBnIj->mBIj')...
+einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3B.PpPhHh,'mneF,FeBnIj->mBIj')...
-0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3B.PPPhHh,'mnEF,EFBnIj->mBIj')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3B.ppPHHh,'mNef,efBINj->mBIj')...
-einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3B.PpPHHh,'mNeF,FeBINj->mBIj')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3B.PPPHHh,'mNEF,EFBINj->mBIj');

Vt3B.oPHh = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(hA,hB,pA,:),T3B.PpPhHh,'nmfe,AfBnIm->ABIe')...
-einsum_kg(sys.vB_oovv(hA,hB,PA,:),T3B.PPPhHh,'nmFe,AFBnIm->ABIe')...
+einsum_kg(sys.vB_oovv(HA,hB,pA,:),T3B.PpPHHh,'Nmfe,AfBINm->ABIe')...
+einsum_kg(sys.vB_oovv(HA,hB,PA,:),T3B.PPPHHh,'NmFe,AFBINm->ABIe')...
-einsum_kg(sys.vB_oovv(hA,HB,pA,:),T3B.PpPhHH,'nMfe,AfBnIM->ABIe')...
-einsum_kg(sys.vB_oovv(hA,HB,PA,:),T3B.PPPhHH,'nMFe,AFBnIM->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,pA,:),T3B.PpPHHH,'NMfe,AfBINM->ABIe')...
+einsum_kg(sys.vB_oovv(HA,HB,PA,:),T3B.PPPHHH,'NMFe,AFBINM->ABIe');

d2 = -0.5*einsum_kg(sys.vC_oovv(hB,hB,pB,:),T3C.PPpHhh,'nmfe,ABfInm->ABIe')...
+0.5*einsum_kg(sys.vC_oovv(hB,hB,PB,:),T3C.PPPHhh,'nmFe,AFBInm->ABIe')...
-einsum_kg(sys.vC_oovv(hB,HB,pB,:),T3C.PPpHhH,'nMfe,ABfInM->ABIe')...
+einsum_kg(sys.vC_oovv(hB,HB,PB,:),T3C.PPPHhH,'nMFe,AFBInM->ABIe')...
-0.5*einsum_kg(sys.vC_oovv(HB,HB,pB,:),T3C.PPpHHH,'NMfe,ABfINM->ABIe')...
+0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),T3C.PPPHHH,'NMFe,AFBINM->ABIe');

Vt3B.PPHv = -d1 - d2;

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

d1 = -einsum_kg(sys.vB_oovv(hA,:,pA,pB),T3B.ppphHh,'nmfe,afenIj->amIj')...
-einsum_kg(sys.vB_oovv(hA,:,pA,PB),T3B.ppPhHh,'nmfE,afEnIj->amIj')...
+einsum_kg(sys.vB_oovv(hA,:,PA,pB),T3B.PpphHh,'nmFe,FaenIj->amIj')...
+einsum_kg(sys.vB_oovv(hA,:,PA,PB),T3B.PpPhHh,'nmFE,FaEnIj->amIj')...
+einsum_kg(sys.vB_oovv(HA,:,pA,pB),T3B.pppHHh,'Nmfe,afeINj->amIj')...
+einsum_kg(sys.vB_oovv(HA,:,pA,PB),T3B.ppPHHh,'NmfE,afEINj->amIj')...
-einsum_kg(sys.vB_oovv(HA,:,PA,pB),T3B.PppHHh,'NmFe,FaeINj->amIj')...
-einsum_kg(sys.vB_oovv(HA,:,PA,PB),T3B.PpPHHh,'NmFE,FaEINj->amIj');

d2 = +0.5*einsum_kg(sys.vC_oovv(hB,:,pB,pB),T3C.pppHhh,'nmfe,afeInj->amIj')...
+einsum_kg(sys.vC_oovv(hB,:,PB,pB),T3C.pPpHhh,'nmFe,aFeInj->amIj')...
+0.5*einsum_kg(sys.vC_oovv(hB,:,PB,PB),T3C.pPPHhh,'nmFE,aFEInj->amIj')...
-0.5*einsum_kg(sys.vC_oovv(HB,:,pB,pB),T3C.pppHhH,'Nmfe,afeIjN->amIj')...
-einsum_kg(sys.vC_oovv(HB,:,PB,pB),T3C.pPpHhH,'NmFe,aFeIjN->amIj')...
-0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),T3C.pPPHhH,'NmFE,aFEIjN->amIj');

Vt3B.poHh = d1 + d2;

d1 = -einsum_kg(sys.vB_oovv(hA,hB,:,pB),T3B.PpphHh,'mnef,AbfmIn->AbIe')...
-einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3B.PpPhHh,'mneF,AbFmIn->AbIe')...
-einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3B.PpphHH,'mNef,AbfmIN->AbIe')...
-einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3B.PpPhHH,'mNeF,AbFmIN->AbIe')...
+einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3B.PppHHh,'Mnef,AbfIMn->AbIe')...
+einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3B.PpPHHh,'MneF,AbFIMn->AbIe')...
+einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3B.PppHHH,'MNef,AbfIMN->AbIe')...
+einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3B.PpPHHH,'MNeF,AbFIMN->AbIe');

d2 = +0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3A.PpphhH,'mnef,AbfmnI->AbIe')...
-0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3A.PPphhH,'mneF,AFbmnI->AbIe')...
+einsum_kg(sys.vA_oovv(HA,hA,:,pA),T3A.PpphHH,'Mnef,AbfnIM->AbIe')...
-einsum_kg(sys.vA_oovv(HA,hA,:,PA),T3A.PPphHH,'MneF,AFbnIM->AbIe')...
+0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3A.PppHHH,'MNef,AbfIMN->AbIe')...
-0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3A.PPpHHH,'MNeF,AFbIMN->AbIe');

Vt3A.PpHv = -d1 - d2;

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
              

d1 = +einsum_kg(sys.vB_oovv(:,hB,pA,pB),T3B.pppHHh,'mnef,aefIJn->amIJ')...
-einsum_kg(sys.vB_oovv(:,hB,PA,pB),T3B.PppHHh,'mnEf,EafIJn->amIJ')...
+einsum_kg(sys.vB_oovv(:,hB,pA,PB),T3B.ppPHHh,'mneF,aeFIJn->amIJ')...
-einsum_kg(sys.vB_oovv(:,hB,PA,PB),T3B.PpPHHh,'mnEF,EaFIJn->amIJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,pB),T3B.pppHHH,'mNef,aefIJN->amIJ')...
-einsum_kg(sys.vB_oovv(:,HB,PA,pB),T3B.PppHHH,'mNEf,EafIJN->amIJ')...
+einsum_kg(sys.vB_oovv(:,HB,pA,PB),T3B.ppPHHH,'mNeF,aeFIJN->amIJ')...
-einsum_kg(sys.vB_oovv(:,HB,PA,PB),T3B.PpPHHH,'mNEF,EaFIJN->amIJ');

d2 = +0.5*einsum_kg(sys.vA_oovv(:,hA,pA,pA),T3A.ppphHH,'mnef,aefnIJ->amIJ')...
+einsum_kg(sys.vA_oovv(:,hA,pA,PA),T3A.PpphHH,'mneF,FaenIJ->amIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,hA,PA,PA),T3A.PPphHH,'mnEF,EFanIJ->amIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,pA,pA),T3A.pppHHH,'mNef,aefIJN->amIJ')...
+einsum_kg(sys.vA_oovv(:,HA,pA,PA),T3A.PppHHH,'mNeF,FaeIJN->amIJ')...
+0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),T3A.PPpHHH,'mNEF,EFaIJN->amIJ');

Vt3A.poHH = d1 + d2;
end