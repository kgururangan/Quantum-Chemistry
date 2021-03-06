function [t3b] = update_t3b_ccsdt3(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, sys, shift)
    
    % get CCSD HBar intermediates 
    [H1A, H1B, H2A, H2B, H2C, VTA, VTB] = get_t3b_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);
    
    PA = sys.PA; PB = sys.PB; HA = sys.HA; HB = sys.HB;
    
    % MM23B + (V*T2*T3)_C    
    I2B_ovoo = H2B.ovoo - einsum_kg(H1A.ov,t2b,'me,ecjk->mcjk');
    I2B_vooo = H2B.vooo - einsum_kg(H1B.ov,t2b,'me,aeik->amik');
    I2A_vooo = H2A.vooo - einsum_kg(H1A.ov,t2a,'me,aeij->amij');    
   
    M23_D1 =  einsum_kg(H2B.vvvo(PA,PB,:,HB) + VTB.vvvo,t2a(PA,:,HA,HA),'BCeK,AeIJ->ABCIJK');
    M23_D2 = -einsum_kg(I2B_ovoo(:,PB,HA,HB) + VTB.ovoo,t2a(PA,PA,HA,:),'mCJK,ABIm->ABCIJK');
    M23_D3 = +einsum_kg(H2B.vvov(PA,PB,HA,:) + VTB.vvov,t2b(PA,:,HA,HB),'ACIe,BeJK->ABCIJK');
    M23_D4 = -einsum_kg(I2B_vooo(PA,:,HA,HB) + VTB.vooo,t2b(PA,PB,HA,:),'AmIK,BCJm->ABCIJK');
    M23_D5 = +einsum_kg(H2A.vvov(PA,PA,HA,:) + VTA.vvov,t2b(:,PB,HA,HB),'ABIe,eCJK->ABCIJK');
    M23_D6 = -einsum_kg(I2A_vooo(PA,:,HA,HA) + VTA.vooo,t2b(PA,PB,:,HB),'AmIJ,BCmK->ABCIJK');
    
    M23_34 = M23_D3 + M23_D4;
    M23_34 = M23_34 - permute(M23_34,[2,1,3,4,5,6]) - permute(M23_34,[1,2,3,5,4,6]) + permute(M23_34,[2,1,3,5,4,6]);
    M23_25 = M23_D2 + M23_D5;
    M23_25 = M23_25 - permute(M23_25,[1,2,3,5,4,6]);
    M23_16 = M23_D1 + M23_D6;
    M23_16 = M23_16 - permute(M23_16,[2,1,3,4,5,6]);
        
    MM23B = M23_16 + M23_25 + M23_34;

    % (HBar*T3)_C
    D1 = -einsum_kg(H1A.oo(HA,HA),t3b,'mi,abcmjk->abcijk');
    D2 = -einsum_kg(H1B.oo(HB,HB),t3b,'mk,abcijm->abcijk');
    D3 = einsum_kg(H1A.vv(PA,PA),t3b,'ae,ebcijk->abcijk');
    D4 = einsum_kg(H1B.vv(PB,PB),t3b,'ce,abeijk->abcijk');
    D5 = 0.5*einsum_kg(H2A.oooo(HA,HA,HA,HA),t3b,'mnij,abcmnk->abcijk');
    D6 = einsum_kg(H2B.oooo(HA,HB,HA,HB),t3b,'mnjk,abcimn->abcijk');
    D7 = 0.5*einsum_kg(H2A.vvvv(PA,PA,PA,PA),t3b,'abef,efcijk->abcijk');
    D8 = einsum_kg(H2B.vvvv(PA,PB,PA,PB),t3b,'bcef,aefijk->abcijk');
    D9 = einsum_kg(H2A.voov(PA,HA,HA,PA),t3b,'amie,ebcmjk->abcijk');    
    D10 = einsum_kg(H2B.voov(PA,HB,HA,PB),t3c,'amie,becjmk->abcijk');    
    D11 = einsum_kg(H2B.ovvo(HA,PB,PA,HB),t3a,'mcek,abeijm->abcijk');
    D12 = einsum_kg(H2C.voov(PB,HB,HB,PB),t3b,'cmke,abeijm->abcijk');
    D13 = -einsum_kg(H2B.vovo(PA,HB,PA,HB),t3b,'amek,ebcijm->abcijk');
    D14 = -einsum_kg(H2B.ovov(HA,PB,HA,PB),t3b,'mcie,abemjk->abcijk');
    
    % diagrams that have A(ab)A(ij)
    D_9_10 = D9 + D10;
    D_9_10 = D_9_10 - permute(D_9_10,[2,1,3,4,5,6]) - permute(D_9_10,[1,2,3,5,4,6]) + permute(D_9_10,[2,1,3,5,4,6]);
    % diagrams that have A(ab)
    D_3_8_13 = D3 + D8 + D13;
    D_3_8_13 = D_3_8_13 - permute(D_3_8_13,[2,1,3,4,5,6]);
    % diagrams that have A(ij)
    D_1_6_14 = D1 + D6 + D14;
    D_1_6_14 = D_1_6_14 - permute(D_1_6_14,[1,2,3,5,4,6]);
     
    X3B_abcijk = MM23B + D2 + D4 + D5 + D7 + D11 + D12 + D_3_8_13 + D_1_6_14 + D_9_10;


     for a = 1:sys.Nact_p_alpha
        for b = a+1:sys.Nact_p_alpha
            for c = 1:sys.Nact_p_beta
                for i = 1:sys.Nact_h_alpha
                    for j = i+1:sys.Nact_h_alpha
                        for k = 1:sys.Nact_h_beta

                            t3b(a,b,c,i,j,k) = t3b(a,b,c,i,j,k) + X3B_abcijk(a,b,c,i,j,k)/...
                                (sys.fa_HH(i,i)+sys.fa_HH(j,j)+sys.fb_HH(k,k)-sys.fa_PP(a,a)-sys.fa_PP(b,b)-sys.fb_PP(c,c)-shift);
                            t3b(b,a,c,i,j,k) = -t3b(a,b,c,i,j,k);
                            t3b(a,b,c,j,i,k) = -t3b(a,b,c,i,j,k);
                            t3b(b,a,c,j,i,k) = t3b(a,b,c,i,j,k);

                        end
                    end
                end
            end
        end
    end


end

function [H1A,H1B,H2A,H2B,H2C,VTA,VTB] = get_t3b_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys)

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

    % (VT3)_C intermediates
    PA = sys.PA; PB  = sys.PB; HA = sys.HA; HB = sys.HB;
        
    VTA.vvov = -0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),t3a,'mnef,abfimn->abie')...
                 -einsum_kg(sys.vB_oovv(HA,HB,:,PB),t3b,'mnef,abfimn->abie');
             
    VTA.vooo = 0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),t3a,'mnef,aefijn->amij')...
                 +einsum_kg(sys.vB_oovv(:,HB,PA,PB),t3b,'mnef,aefijn->amij');
             
    VTB.vvvo = -0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),t3b,'mnef,afbmnj->abej')...
                 -einsum_kg(sys.vB_oovv(HA,HB,:,PB),t3c,'mnef,afbmnj->abej');
             
    VTB.ovoo = 0.5*einsum_kg(sys.vA_oovv(:,HA,PA,PA),t3b,'mnef,efbinj->mbij')...
                 +einsum_kg(sys.vB_oovv(:,HB,PA,PB),t3c,'mnef,efbinj->mbij');
             
    VTB.vvov = -einsum_kg(sys.vB_oovv(HA,HB,PA,:),t3b,'nmfe,afbinm->abie')...
                 -0.5*einsum_kg(sys.vC_oovv(HB,HB,PB,:),t3c,'nmfe,afbinm->abie');
             
    VTB.vooo = einsum_kg(sys.vB_oovv(HA,:,PA,PB),t3b,'nmfe,afeinj->amij')...
                 +0.5*einsum_kg(sys.vC_oovv(HB,:,PB,PB),t3c,'nmfe,afeinj->amij');
             
end

