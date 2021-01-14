function [X3D_abcijk,VTC] = build_t3d(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, sys, shift)

    [H1B, H2B, H2C, VTC] = get_t3d_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);

    I2C_vvov = H2C.vvov+einsum_kg(H1B.ov,t2c,'me,abim->abie');
    
    % MM23D
    M23_D1 =   -einsum_kg(H2C.vooo + VTC.vooo,t2c,'amij,bcmk->abcijk'); 
    M23_D1 = M23_D1 - permute(M23_D1,[1,2,3,6,5,4]) - permute(M23_D1,[1,2,3,4,6,5]) - permute(M23_D1,[3,2,1,4,5,6])...
                    - permute(M23_D1,[2,1,3,4,5,6]) + permute(M23_D1,[2,1,3,6,5,4]) + permute(M23_D1,[3,2,1,6,5,4]) ...
                    + permute(M23_D1,[2,1,3,4,6,5]) + permute(M23_D1,[3,2,1,4,6,5]);

         
    M23_D2 =   +einsum_kg(I2C_vvov + VTC.vvov,t2c,'abie,ecjk->abcijk');
    M23_D2 = M23_D2 - permute(M23_D2,[1,2,3,5,4,6]) - permute(M23_D2,[1,2,3,6,5,4]) ...
                    - permute(M23_D2,[3,2,1,4,5,6]) - permute(M23_D2,[1,3,2,4,5,6]) + permute(M23_D2,[3,2,1,5,4,6]) ...
                    + permute(M23_D2,[1,3,2,5,4,6]) + permute(M23_D2,[3,2,1,6,5,4]) + permute(M23_D2,[1,3,2,6,5,4]);
        
    MM23D = M23_D1 + M23_D2;

    % (HBar*T3)_C

    D1 = -einsum_kg(H1B.oo,t3d,'mk,abcijm->abcijk');
    D2 = einsum_kg(H1B.vv,t3d,'ce,abeijk->abcijk');
    D3 = 0.5*einsum_kg(H2C.oooo,t3d,'mnij,abcmnk->abcijk');
    D4 = 0.5*einsum_kg(H2C.vvvv,t3d,'abef,efcijk->abcijk');
    D5 = einsum_kg(H2B.ovvo,t3c,'maei,ebcmjk->abcijk');
    D6 = einsum_kg(H2C.voov,t3d,'amie,ebcmjk->abcijk');

%     % A(k/ij)
%     D13 = D1 + D3;
%     D13 = D13 - permute(D13,[1,2,3,6,5,4]) - permute(D13,[1,2,3,4,6,5]);
%     
%     % A(c/ab)
%     D24 = D2 + D4;
%     D24 = D24 - permute(D24,[3,2,1,4,5,6]) - permute(D24,[1,3,2,4,5,6]);
%    
%     % A(i/jk)A(a/bc)
%     D56 = D5 + D6;
%     D56 = D56 - permute(D56,[1,2,3,5,4,6]) - permute(D56,[1,2,3,6,5,4])...
%                -permute(D56,[2,1,3,4,5,6]) - permute(D56,[3,2,1,4,5,6]) ...
%                +permute(D56,[2,1,3,5,4,6]) + permute(D56,[2,1,3,6,5,4]) ...
%                +permute(D56,[3,2,1,5,4,6]) + permute(D56,[3,2,1,6,5,4]);  

     D1 = D1 - permute(D1,[1,2,3,6,5,4]) - permute(D1,[1,2,3,4,6,5]);
     D2 = D2 - permute(D2,[3,2,1,4,5,6]) - permute(D2,[1,3,2,4,5,6]);
     D3 = D3 - permute(D3,[1,2,3,6,5,4]) - permute(D3,[1,2,3,4,6,5]);
     D4 = D4 - permute(D4,[3,2,1,4,5,6]) - permute(D4,[1,3,2,4,5,6]);
     D5 = D5 - permute(D5,[1,2,3,5,4,6]) - permute(D5,[1,2,3,6,5,4])...
              -permute(D5,[2,1,3,4,5,6]) - permute(D5,[3,2,1,4,5,6]) ...
              +permute(D5,[2,1,3,5,4,6]) + permute(D5,[2,1,3,6,5,4]) ...
              +permute(D5,[3,2,1,5,4,6]) + permute(D5,[3,2,1,6,5,4]); 
     D6 = D6 - permute(D6,[1,2,3,5,4,6]) - permute(D6,[1,2,3,6,5,4])...
              -permute(D6,[2,1,3,4,5,6]) - permute(D6,[3,2,1,4,5,6]) ...
              +permute(D6,[2,1,3,5,4,6]) + permute(D6,[2,1,3,6,5,4]) ...
              +permute(D6,[3,2,1,5,4,6]) + permute(D6,[3,2,1,6,5,4]); 
     
     X3D_abcijk = MM23D + D1 + D2 + D3 + D4 + D5 + D6;
     
       

end

function [H1B, H2B, H2C, VTC] = get_t3d_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys)

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
            
    % VT3 intermediates
    VTC.vvov = -0.5*einsum_kg(sys.vC_oovv,t3d,'mnef,abfimn->abie')...
                 -einsum_kg(sys.vB_oovv,t3c,'nmfe,fabnim->abie');
    VTC.vvvo = -permute(VTC.vvov,[1,2,4,3]);
             
    VTC.vooo = einsum_kg(sys.vB_oovv,t3c,'nmfe,faenij->amij')...
                 +0.5*einsum_kg(sys.vC_oovv,t3d,'mnef,aefijn->amij');
    VTC.ovoo = -permute(VTC.vooo,[2,1,3,4]);
            
end