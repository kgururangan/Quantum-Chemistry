function [X3B,VTA,VTB] = build_t3b(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, sys, shift)
    
    % get CCSD HBar intermediates 
    [H1A, H1B, H2A, H2B, H2C, VTA, VTB] = get_t3b_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);
    
    % MM23B + (V*T2*T3)_C    
    I2B_ovoo = H2B.ovoo - einsum_kg(H1A.ov,t2b,'me,ecjk->mcjk');
    I2B_vooo = H2B.vooo - einsum_kg(H1B.ov,t2b,'me,aeik->amik');
    I2A_vooo = H2A.vooo - einsum_kg(H1A.ov,t2a,'me,aeij->amij');    
%    
    M23_D1 = einsum_kg(H2B.vvvo + VTB.vvvo,t2a,'bcek,aeij->abcijk');
    M23_D2 = -einsum_kg(I2B_ovoo + VTB.ovoo,t2a,'mcjk,abim->abcijk');
    M23_D3 = +einsum_kg(H2B.vvov + VTB.vvov,t2b,'acie,bejk->abcijk');
    M23_D4 = -einsum_kg(I2B_vooo + VTB.vooo,t2b,'amik,bcjm->abcijk');
    M23_D5 = +einsum_kg(H2A.vvov + VTA.vvov,t2b,'abie,ecjk->abcijk');
    M23_D6 = -einsum_kg(I2A_vooo + VTA.vooo,t2b,'amij,bcmk->abcijk');

%     M23_D1 = einsum_kg(H2B.vvvo,t2a,'bcek,aeij->abcijk');
%     M23_D2 = -einsum_kg(I2B_ovoo,t2a,'mcjk,abim->abcijk');
%     M23_D3 = +einsum_kg(H2B.vvov,t2b,'acie,bejk->abcijk');
%     M23_D4 = -einsum_kg(I2B_vooo,t2b,'amik,bcjm->abcijk');
%     M23_D5 = +einsum_kg(H2A.vvov,t2b,'abie,ecjk->abcijk');
%     M23_D6 = -einsum_kg(I2A_vooo,t2b,'amij,bcmk->abcijk');
%     
    M23_34 = M23_D3 + M23_D4;
    M23_34 = M23_34 - permute(M23_34,[2,1,3,4,5,6]) - permute(M23_34,[1,2,3,5,4,6]) + permute(M23_34,[2,1,3,5,4,6]);
    M23_25 = M23_D2 + M23_D5;
    M23_25 = M23_25 - permute(M23_25,[1,2,3,5,4,6]);
    M23_16 = M23_D1 + M23_D6;
    M23_16 = M23_16 - permute(M23_16,[2,1,3,4,5,6]);
        
    MM23B = M23_16 + M23_25 + M23_34;

    % (HBar*T3)_C
    D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
    D2 = -einsum_kg(H1B.oo,t3b,'mk,abcijm->abcijk');
    D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
    D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
    D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
    D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
    D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
    D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk');
    D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');    
    D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk');    
    D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
    D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
    D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
    D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
    
%     % diagrams that have A(ab)A(ij)
%     D_9_10 = D9 + D10;
%     D_9_10 = D_9_10 - permute(D_9_10,[2,1,3,4,5,6]) - permute(D_9_10,[1,2,3,5,4,6]) + permute(D_9_10,[2,1,3,5,4,6]);
%     % diagrams that have A(ab)
%     D_3_8_13 = D3 + D8 + D13;
%     D_3_8_13 = D_3_8_13 - permute(D_3_8_13,[2,1,3,4,5,6]);
%     % diagrams that have A(ij)
%     D_1_6_14 = D1 + D6 + D14;
%     D_1_6_14 = D_1_6_14 - permute(D_1_6_14,[1,2,3,5,4,6]);
     
%     X3B = MM23B + D2 + D4 + D5 + D7 + D11 + D12 + D_3_8_13 + D_1_6_14 + D_9_10;

% projection 4
%     D1 = D1 - permute(D1,[1,2,3,5,4,6]);
%     D3 = D3 - permute(D3,[2,1,3,4,5,6]);
%     D9 = D9 - permute(D9,[2,1,3,4,5,6]) - permute(D9,[1,2,3,5,4,6]) + permute(D9,[2,1,3,5,4,6]);
%     D10 = D10 - permute(D10,[2,1,3,4,5,6]) - permute(D10,[1,2,3,5,4,6]) + permute(D10,[2,1,3,5,4,6]);
%     D13 = D13 - permute(D13,[2,1,3,4,5,6]);
%     D14 = D14 - permute(D14,[1,2,3,5,4,6]);
%     D6 = D6 - permute(D6,[1,2,3,5,4,6]);
%     D8 = D8 - permute(D8,[2,1,3,4,5,6]);
%     
%     M23_D6 = M23_D6 - permute(M23_D6,[2,1,3,4,5,6]);
%     M23_D4 = M23_D4 - permute(M23_D4,[2,1,3,4,5,6]) - permute(M23_D4,[1,2,3,5,4,6]) + permute(M23_D4,[2,1,3,5,4,6]);
%     M23_D2 = M23_D2 - permute(M23_D2,[1,2,3,5,4,6]);
%     M23_D5 = M23_D5 - permute(M23_D5,[1,2,3,5,4,6]);
%     M23_D1 = M23_D1 - permute(M23_D1,[2,1,3,4,5,6]);
%     M23_D3 = M23_D3 - permute(M23_D3,[2,1,3,4,5,6]) - permute(M23_D3,[1,2,3,5,4,6]) + permute(M23_D3,[2,1,3,5,4,6]);
%     
%     X3B = D1+D2+D3+D4+D9+D11+D12+D10+D13+D14+D5+D6+D7+D8...
%           +M23_D6+M23_D4+M23_D2+M23_D5+M23_D1+M23_D3;


% projection check
    D1 = D1 - permute(D1,[1,2,3,5,4,6]);
    D3 = D3 - permute(D3,[2,1,3,4,5,6]);
    D9 = D9 - permute(D9,[2,1,3,4,5,6]) - permute(D9,[1,2,3,5,4,6]) + permute(D9,[2,1,3,5,4,6]);
    D10 = D10 - permute(D10,[2,1,3,4,5,6]) - permute(D10,[1,2,3,5,4,6]) + permute(D10,[2,1,3,5,4,6]);
    D13 = D13 - permute(D13,[2,1,3,4,5,6]);
    D14 = D14 - permute(D14,[1,2,3,5,4,6]);
    D6 = D6 - permute(D6,[1,2,3,5,4,6]);
    D8 = D8 - permute(D8,[2,1,3,4,5,6]);

    X3B = MM23B + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12 + D13 + D14;
    %X3B = D1;





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
    VTA.vvov = -0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,abfimn->abie')...
                 -einsum_kg(sys.vB_oovv,t3b,'mnef,abfimn->abie');
             
    VTA.vooo = 0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,aefijn->amij')...
                 +einsum_kg(sys.vB_oovv,t3b,'mnef,aefijn->amij');
             
    VTB.vvvo = -0.5*einsum_kg(sys.vA_oovv,t3b,'mnef,afbmnj->abej')...
                 -einsum_kg(sys.vB_oovv,t3c,'mnef,afbmnj->abej');
             
    VTB.ovoo = 0.5*einsum_kg(sys.vA_oovv,t3b,'mnef,efbinj->mbij')...
                 +einsum_kg(sys.vB_oovv,t3c,'mnef,efbinj->mbij');
             
    VTB.vvov = -einsum_kg(sys.vB_oovv,t3b,'nmfe,afbinm->abie')...
                 -0.5*einsum_kg(sys.vC_oovv,t3c,'nmfe,afbinm->abie');
             
    VTB.vooo = einsum_kg(sys.vB_oovv,t3b,'nmfe,afeinj->amij')...
                 +0.5*einsum_kg(sys.vC_oovv,t3c,'nmfe,afeinj->amij');
             
end

% DEBUG BUILD
%   %
%     M23_D1 = einsum_kg(H2B.vvvo,t2a,'bcek,aeij->abcijk');
% 
%     M23_D2 = -einsum_kg(I2B_ovoo,t2a,'mcjk,abim->abcijk');
%     
%     M23_D3 = +einsum_kg(H2B.vvov,t2b,'acie,bejk->abcijk');
% 
%     M23_D4 = -einsum_kg(I2B_vooo,t2b,'amik,bcjm->abcijk');
% 
%     M23_D5 = +einsum_kg(H2A.vvov,t2b,'abie,ecjk->abcijk');
% 
%     M23_D6 = -einsum_kg(I2A_vooo,t2b,'amij,bcmk->abcijk');
%
%     M23_34 = M23_D3 + M23_D4;
%     M23_34 = M23_34 - permute(M23_34,[2,1,3,4,5,6]) - permute(M23_34,[1,2,3,5,4,6]) + permute(M23_34,[2,1,3,5,4,6]);
%         
%     M23_25 = M23_D2 + M23_D5;
%     M23_25 = M23_25 - permute(M23_25,[1,2,3,5,4,6]);
% 
%     M23_16 = M23_D1 + M23_D6;
%     M23_16 = M23_16 - permute(M23_16,[2,1,3,4,5,6]);
%         
%     MM23B = M23_16 + M23_25 + M23_34;
%    
%     VTA_vvov = -0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,abfimn->abie') - einsum_kg(sys.vB_oovv,t3b,'mnef,abfimn->abie');                     
%     VT_1 = einsum_kg(VTA_vvov,t2b,'abie,ecjk->abcijk');
%     VT_1 = VT_1 - permute(VT_1,[1,2,3,5,4,6]);
%              
%     VTB_vvvo = -0.5*einsum_kg(sys.vA_oovv,t3b,'mnef,bfcmnk->bcek') - einsum_kg(sys.vB_oovv,t3c,'mnef,bfcmnk->bcek');
%     VT_2 = einsum_kg(VTB_vvvo,t2a,'bcek,aeij->abcijk');
%     VT_2 = VT_2 - permute(VT_2,[2,1,3,4,5,6]);
%     
%     tmp = -einsum_kg(sys.vB_oovv,t3b,'nmfe,bfcjnm->bcje') - 0.5*einsum_kg(sys.vC_oovv,t3c,'nmfe,bfcjnm->bcje');
%     VT_3 = einsum_kg(tmp,t2b,'bcje,aeik->abcijk');
%     VT_3 = VT_3 - permute(VT_3,[2,1,3,4,5,6]) - permute(VT_3,[1,2,3,5,4,6]) + permute(VT_3,[2,1,3,5,4,6]);
%         
%     tmp = 0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,ebfijn->mbij') + einsum_kg(sys.vB_oovv,t3b,'mnef,ebfijn->mbij');
%     VT_4 = -einsum_kg(tmp,t2b,'mbij,acmk->abcijk');
%     VT_4 = VT_4 - permute(VT_4,[2,1,3,4,5,6]);
%     
%     tmp = 0.5*einsum_kg(sys.vA_oovv,t3b,'mnef,efcink->mcik') + einsum_kg(sys.vB_oovv,t3c,'mnef,efcink->mcik');
%     VT_5 = -einsum_kg(tmp,t2a,'mcik,abmj->abcijk');
%     VT_5 = VT_5 - permute(VT_5,[1,2,3,5,4,6]);
%     
%     tmp = einsum_kg(sys.vB_oovv,t3b,'nmfe,afeink->amik') + 0.5*einsum_kg(sys.vC_oovv,t3c,'nmfe,afeink->amik');
%     VT_6 = -einsum_kg(tmp,t2b,'amik,bcjm->abcijk');
%     VT_6 = VT_6 - permute(VT_6,[2,1,3,4,5,6]) - permute(VT_6,[1,2,3,5,4,6]) + permute(VT_6,[2,1,3,5,4,6]);
%     
%     MM23B = MM23B + VT_1 + VT_2 + VT_3 + VT_4 + VT_5 + VT_6;
%
%
%     D1 = -einsum_kg(H1A.oo-diag(diag(sys.fa_oo)),t3b,'mi,abcmjk->abcijk');
%     D1 = D1 - permute(D1,[1,2,3,5,4,6]);
%     
%     D2 = -einsum_kg(H1B.oo-diag(diag(sys.fb_oo)),t3b,'mk,abcijm->abcijk');
%     
%     D3 = einsum_kg(H1A.vv-diag(diag(sys.fa_vv)),t3b,'ae,ebcijk->abcijk');
%     D3 = D3 - permute(D3,[2,1,3,4,5,6]);
%     
%     D4 = einsum_kg(H1B.vv-diag(diag(sys.fb_vv)),t3b,'ce,abeijk->abcijk');
%     
%     D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
%     
%     D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
%     D6 = D6 - permute(D6,[1,2,3,5,4,6]);
%     
%     D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
%     
%     D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk');
%     D8 = D8 - permute(D8,[2,1,3,4,5,6]);
%     
%     D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');
%     D9 = D9 - permute(D9,[1,2,3,5,4,6]) - permute(D9,[2,1,3,4,5,6]) + permute(D9,[2,1,3,5,4,6]);
%     
%     D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk');
%     D10 = D10 - permute(D10,[1,2,3,5,4,6]) - permute(D10,[2,1,3,4,5,6]) + permute(D10,[2,1,3,5,4,6]);
%     
%     D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
%     
%     D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
%     
%     D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
%     D13 = D13 - permute(D13,[2,1,3,4,5,6]);
%     
%     D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
%     D14 = D14 - permute(D14,[1,2,3,5,4,6]);
%     
%     X3B_abcijk = MM23B + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 ...
%                        + D11 + D12 + D13 + D14;
