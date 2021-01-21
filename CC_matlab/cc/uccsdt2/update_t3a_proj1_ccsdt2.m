function [t3a] = update_t3a_proj1_ccsdt2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift)

    H1A = HBar_t.H1A; H2A = HBar_t.H2A; H2B = HBar_t.H2B;
    [Vt3A] = get_t3a_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys);
    
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

    I2A_vvov = H2A.vvov+einsum_kg(H1A.ov,t2a,'me,abim->abie');
    
    % MM23A + (V*T2*T3)_C
    M23_D1 =   -einsum_kg(H2A.vooo(PA,:,HA,HA) + Vt3A.PoHH,t2a(PA,PA,:,HA),'AmIJ,BCmK->ABCIJK'); 
    M23_D1 = M23_D1 - permute(M23_D1,[1,2,3,6,5,4]) - permute(M23_D1,[1,2,3,4,6,5]) - permute(M23_D1,[3,2,1,4,5,6])...
                    - permute(M23_D1,[2,1,3,4,5,6]) + permute(M23_D1,[2,1,3,6,5,4]) + permute(M23_D1,[3,2,1,6,5,4]) ...
                    + permute(M23_D1,[2,1,3,4,6,5]) + permute(M23_D1,[3,2,1,4,6,5]);

         
    M23_D2 =   +einsum_kg(I2A_vvov(PA,PA,HA,:) + Vt3A.PPHv,t2a(:,PA,HA,HA),'ABIe,eCJK->ABCIJK');
    M23_D2 = M23_D2 - permute(M23_D2,[1,2,3,5,4,6]) - permute(M23_D2,[1,2,3,6,5,4]) ...
                    - permute(M23_D2,[3,2,1,4,5,6]) - permute(M23_D2,[1,3,2,4,5,6]) + permute(M23_D2,[3,2,1,5,4,6]) ...
                    + permute(M23_D2,[1,3,2,5,4,6]) + permute(M23_D2,[3,2,1,6,5,4]) + permute(M23_D2,[1,3,2,6,5,4]);
                
    MM23A = M23_D1 + M23_D2;

    % (HBar*T3)_C    
    %D1 = -einsum_kg(H1A.oo(HA,HA),T3A.PPPHHH,'mk,abcijm->abcijk');
 D1 = -einsum_kg(H1A.oo(hA,HA),T3A.PPPhHH,'mK,ABCmIJ->ABCIJK')...
-einsum_kg(H1A.oo(HA,HA),T3A.PPPHHH,'MK,ABCIJM->ABCIJK');   

%D2 = einsum_kg(H1A.vv(PA,PA),T3A.PPPHHH,'ce,abeijk->abcijk');
D2 = +einsum_kg(H1A.vv(PA,pA),T3A.PPpHHH,'Ce,ABeIJK->ABCIJK')...
+einsum_kg(H1A.vv(PA,PA),T3A.PPPHHH,'CE,ABEIJK->ABCIJK');

%D3 = 0.5*einsum_kg(H2A.oooo(HA,HA,HA,HA),T3A.PPPHHH,'mnij,abcmnk->abcijk'); 
D3 = +0.5*einsum_kg(H2A.oooo(hA,hA,HA,HA),T3A.PPPhhH,'mnIJ,ABCmnK->ABCIJK')...
-einsum_kg(H2A.oooo(HA,hA,HA,HA),T3A.PPPhHH,'MnIJ,ABCnMK->ABCIJK')...
+0.5*einsum_kg(H2A.oooo(HA,HA,HA,HA),T3A.PPPHHH,'MNIJ,ABCMNK->ABCIJK');
  
%D4 = 0.5*einsum_kg(H2A.vvvv(PA,PA,PA,PA),T3A.PPPHHH,'abef,efcijk->abcijk');
D4 = +0.5*einsum_kg(H2A.vvvv(PA,PA,pA,pA),T3A.PppHHH,'ABef,CefIJK->ABCIJK')...
-einsum_kg(H2A.vvvv(PA,PA,PA,pA),T3A.PPpHHH,'ABEf,ECfIJK->ABCIJK')...
+0.5*einsum_kg(H2A.vvvv(PA,PA,PA,PA),T3A.PPPHHH,'ABEF,EFCIJK->ABCIJK');

% D5 = einsum_kg(H2A.voov(PA,HA,HA,PA),T3A.PPPHHH,'cmke,abeijm->abcijk');
D5 = +einsum_kg(H2A.voov(PA,hA,HA,pA),T3A.PPphHH,'CmKe,ABemIJ->ABCIJK')...
+einsum_kg(H2A.voov(PA,hA,HA,PA),T3A.PPPhHH,'CmKE,ABEmIJ->ABCIJK')...
+einsum_kg(H2A.voov(PA,HA,HA,pA),T3A.PPpHHH,'CMKe,ABeIJM->ABCIJK')...
+einsum_kg(H2A.voov(PA,HA,HA,PA),T3A.PPPHHH,'CMKE,ABEIJM->ABCIJK');
  
%D6 = einsum_kg(H2B.voov(PA,HB,HA,PB),T3B.PPPHHH,'cmke,abeijm->abcijk');
D6 = +einsum_kg(H2B.voov(PA,hB,HA,pB),T3B.PPpHHh,'CmKe,ABeIJm->ABCIJK')...
+einsum_kg(H2B.voov(PA,hB,HA,PB),T3B.PPPHHh,'CmKE,ABEIJm->ABCIJK')...
+einsum_kg(H2B.voov(PA,HB,HA,pB),T3B.PPpHHH,'CMKe,ABeIJM->ABCIJK')...
+einsum_kg(H2B.voov(PA,HB,HA,PB),T3B.PPPHHH,'CMKE,ABEIJM->ABCIJK');
    
    % A(k/ij)
    D13 = D1 + D3;
    D13 = D13 - permute(D13,[1,2,3,6,5,4]) - permute(D13,[1,2,3,4,6,5]);
    
    % A(c/ab)
    D24 = D2 + D4;
    D24 = D24 - permute(D24,[3,2,1,4,5,6]) - permute(D24,[1,3,2,4,5,6]);
   
    % A(k/ij)A(c/ab)
    D56 = D5 + D6;
    D56 = D56 - permute(D56,[1,2,3,6,5,4]) - permute(D56,[1,2,3,4,6,5])...
              - permute(D56,[3,2,1,4,5,6]) - permute(D56,[1,3,2,4,5,6]) ...
              + permute(D56,[3,2,1,6,5,4]) + permute(D56,[3,2,1,4,6,5]) ...
              + permute(D56,[1,3,2,6,5,4]) + permute(D56,[1,3,2,4,6,5]);    
     
     X3A_PPPHHH = MM23A + D13 + D24 + D56;
     
     t3a = T3A.PPPHHH;
     for a = 1:sys.Nact_p_alpha
        for b = a+1:sys.Nact_p_alpha
            for c = b+1:sys.Nact_p_alpha
                for i = 1:sys.Nact_h_alpha
                    for j = i+1:sys.Nact_h_alpha
                        for k = j+1:sys.Nact_h_alpha
                            
                            % (1)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(a,b,c,i,j,k) = t3a(a,b,c,i,j,k) + X3A_PPPHHH(a,b,c,i,j,k)/...
                                (sys.fa_HH(i,i)+sys.fa_HH(j,j)+sys.fa_HH(k,k)-sys.fa_PP(a,a)-sys.fa_PP(b,b)-sys.fa_PP(c,c)-shift);    
                            t3a(a,b,c,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(a,b,c,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(a,b,c,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(a,b,c,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(a,b,c,k,j,i) = -t3a(a,b,c,i,j,k);
                            
                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(b,a,c,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(b,a,c,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(b,a,c,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(a,c,b,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(a,c,b,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(a,c,b,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(c,b,a,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(c,b,a,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(c,b,a,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(b,c,a,i,j,k) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(b,c,a,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(b,c,a,k,j,i) = -t3a(a,b,c,i,j,k);
                            
                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(c,a,b,i,j,k) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(c,a,b,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(c,a,b,k,j,i) = -t3a(a,b,c,i,j,k);
                        end
                    end
                end
            end
        end
    end
    %T3A.PPPHHH = t3a;   

end

function [Vt3A] = get_t3a_intermediates(t1a,t1b,t2a,t2b,t2c,T3A,T3B,T3C,T3D,sys)

%     % h1A(me)
%     H1A.ov = sys.fa_ov + einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me') + einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');
% 
%     % h1A(mi)
%     H1A.oo = sys.fa_oo_masked ...
%              +einsum_kg(sys.fa_ov,t1a,'me,ei->mi')...
%              +einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi')...
%              +einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi')...
%              +0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi')...
%              +einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi')...
%              +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t1a,'me,ei->mi')...
%              +einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t1a,'me,ei->mi');
%          
%     % h1A(ae)
%     H1A.vv = sys.fa_vv_masked ...
%              -einsum_kg(sys.fa_ov,t1a,'me,am->ae')...
%              +einsum_kg(sys.vA_vovv,t1a,'amef,fm->ae')...
%              +einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae')...
%              -0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,afmn->ae')...
%              -einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ae')...
%              -einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me'),t1a,'me,am->ae')...
%              -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'),t1a,'me,am->ae');
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
%     % h2B(amie)
%     H2B.voov = sys.vB_voov ...
%                -einsum_kg(sys.vB_ooov,t1a,'nmie,an->amie')...
%                +einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie')...
%                -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie')...
%                +einsum_kg(sys.vB_oovv,t2a,'nmfe,afin->amie')...
%                +einsum_kg(sys.vC_oovv,t2b,'mnef,afin->amie');
           
    % (VT3)_C 
    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
        
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
end