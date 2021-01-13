function [t2c] = update_t2c(t1a, t1b, t2a, t2b, t2c, sys, shift)
   
    % H(CCS) intermediates          
    h1B_ov = sys.fb_ov...
             +einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me')...
             +einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
         
    h1B_oo = sys.fb_oo...
             +einsum_kg(h1B_ov,t1b,'me,ei->mi')...
             +einsum_kg(sys.vB_oovo,t1a,'nmfi,fn->mi')...
             +einsum_kg(sys.vC_ooov,t1b,'mnif,fn->mi');
         
    h1B_vv = sys.fb_vv...
             -einsum_kg(h1B_ov,t1b,'me,am->ae')...
             +einsum_kg(sys.vB_ovvv,t1a,'nafe,fn->ae')...
             +einsum_kg(sys.vC_vovv,t1b,'anef,fn->ae');
          
    h2C_oooo = sys.vC_oooo...
                +einsum_kg(sys.vC_ooov,t1b,'mnie,ej->mnij')...
	    			-einsum_kg(sys.vC_ooov,t1b,'mnje,ei->mnij')...
                +einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ei->mnif'),t1b,'mnif,fj->mnij'); 
            
    h2C_vvvv = sys.vC_vvvv...
                -einsum_kg(sys.vC_ovvv,t1b,'mbef,am->abef')...
	     			+einsum_kg(sys.vC_ovvv,t1b,'maef,bm->abef')...
                +einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,bn->mbef'),t1b,'mbef,am->abef');
                
    h2B_ovvo = sys.vB_ovvo ...
               +einsum_kg(sys.vB_ovvv,t1b,'maef,fi->maei')...
               -einsum_kg(sys.vB_oovo,t1b,'mnei,an->maei')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->maei');
           
    h2C_voov = sys.vC_voov ...
               -einsum_kg(sys.vC_oovo,t1b,'mnei,an->amie')...
               +einsum_kg(sys.vC_vovv,t1b,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fi->mnei'),t1b,'mnei,an->amie');
           
    % <ijab|H(CCS)|0> intermediates - need weights of 0.5 to compensate diagrams that do not have A(ij) or A(ab)
    I2C_ovoo = sys.vC_ovoo ...
               -0.5*einsum_kg(sys.vC_oooo,t1b,'mnij,bn->mbij')...
               +einsum_kg(einsum_kg(sys.vC_ovvv,t1b,'mbef,ei->mbif'),t1b,'mbif,fj->mbij')...
               -0.5*einsum_kg(einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,fj->mnej'),t1b,'mnej,ei->mnij'),t1b,'mnij,bn->mbij')...
               +einsum_kg(sys.vC_ovov,t1b,'mbif,fj->mbij')...
                    -einsum_kg(sys.vC_ovov,t1b,'mbjf,fi->mbij');
                
    I2C_vvvo = sys.vC_vvvo...
               +0.5*einsum_kg(sys.vC_vvvv,t1b,'abef,fj->abej')...
               +einsum_kg(einsum_kg(sys.vC_oovo,t1b,'mnej,am->anej'),t1b,'anej,bn->abej');
             
    % VT2 intermediates - need weights of 0.5 to compensate diagrams that do not have A(ij) or A(ab)
    VT2_1B_oo = 0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,efin->mi') + einsum_kg(sys.vB_oovv,t2b,'nmfe,feni->mi');
    VT2_1B_vv = -0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,afmn->ae') - einsum_kg(sys.vB_oovv,t2b,'nmfe,fanm->ae');
    VT2_2C_oooo = 0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,efij->mnij');
    VT2_2B_ovvo = einsum_kg(sys.vB_oovv,t2c,'mnef,afin->maei') + 0.5*einsum_kg(sys.vA_oovv,t2b,'mnef,fani->maei');
    VT2_2C_voov = 0.5*einsum_kg(sys.vC_oovv,t2c,'mnef,afin->amie');
    
    % contraction intermdiates
    I1B_oo = h1B_oo + VT2_1B_oo;
    I1B_vv = h1B_vv + VT2_1B_vv;
    I2C_voov = h2C_voov + VT2_2C_voov;
    I2B_ovvo = h2B_ovvo + VT2_2B_ovvo;
    I2C_oooo = h2C_oooo + VT2_2C_oooo;
    
    % <ijab| H(CCS) | 0> + <ijab| (H(CCS)T2)_C |0> + <ijab| 0.5*(H(CCS)T2^2)_C |0>
    D1 = -einsum_kg(I2C_ovoo,t1b,'mbij,am->abij'); 
    
    D2 = einsum_kg(I2C_vvvo,t1b,'abej,ei->abij');
    
    D3 = einsum_kg(I1B_vv,t2c,'ae,ebij->abij');
    
    D4 = -einsum_kg(I1B_oo,t2c,'mi,abmj->abij');
    
    D5 = einsum_kg(I2C_voov,t2c,'amie,ebmj->abij');
    
    D6 = einsum_kg(I2B_ovvo,t2b,'maei,ebmj->abij');
    
    D7 = 0.5*einsum_kg(h2C_vvvv,t2c,'abef,efij->abij');
    
    D8 = 0.5*einsum_kg(I2C_oooo,t2c,'mnij,abmn->abij');
    
    % diagrams that have A(ab)
    D13 = D1 + D3;
    D13 = D13 - permute(D13,[2,1,3,4]);
    
    % diagrams that have A(ij)
    D24 = D2 + D4;
    D24 = D24 - permute(D24,[1,2,4,3]);
        
    % diagrams that have A(ab)A(ij)
    D56 = D5 + D6;
    D56 = D56 - permute(D56,[2,1,3,4]) - permute(D56,[1,2,4,3]) + permute(D56,[2,1,4,3]);
    
    % total contribution
    X2C = sys.vC_vvoo + D13 + D24 + D56 + D7 + D8;


    for i = 1:sys.Nocc_beta
        for j = i+1:sys.Nocc_beta
            for a = 1:sys.Nvir_beta
                for b = a+1:sys.Nvir_beta
                    t2c(a,b,i,j) = t2c(a,b,i,j) + X2C(a,b,i,j)/(sys.fb_oo(i,i)+sys.fb_oo(j,j)-sys.fb_vv(a,a)-sys.fb_vv(b,b)-shift);                
                    t2c(b,a,i,j) = -t2c(a,b,i,j);
                    t2c(a,b,j,i) = -t2c(a,b,i,j);
                    t2c(b,a,j,i) = t2c(a,b,i,j);
                end
            end
        end
    end


end



%     d1 = sys.vC_vvoo;
%     d2 = einsum_kg(sys.vC_vvvo,t1b,'abej,ei->abij');
%     d2 = d2 - permute(d2,[1,2,4,3]);
%     d3 = -einsum_kg(sys.vC_ovoo,t1b,'mbij,am->abij');
%     d3 = d3 - permute(d3,[2,1,3,4]);
%     d4 = einsum_kg(einsum_kg(sys.vC_vvvv,t1b,'abef,fj->abej'),t1b,'abej,ei->abij');
%     d5 = einsum_kg(einsum_kg(sys.vC_oooo,t1b,'mnij,am->anij'),t1b,'anij,bn->abij');
%     d6 = einsum_kg(einsum_kg(einsum_kg(sys.vC_ooov,t1b,'mnie,am->anie'),t1b,'anie,bn->abie'),t1b,'abie,ej->abij');
%     d6 = d6 - permute(d6,[1,2,4,3]);
%     d7 = -einsum_kg(einsum_kg(einsum_kg(sys.vC_ovvv,t1b,'mbef,ei->mbif'),t1b,'mbif,fj->mbij'),t1b,'mbij,am->abij');
%     d7 = d7 - permute(d7,[2,1,3,4]);
%     d8 = einsum_kg(einsum_kg(einsum_kg(einsum_kg(sys.vC_oovv,t1b,'mnef,ei->mnif'),t1b,'mnif,fj->mnij'),t1b,'mnij,bn->mbij'),t1b,'mbij,am->abij');
%     d9 = -einsum_kg(einsum_kg(sys.vC_ovvo,t1b,'mbej,ei->mbij'),t1b,'mbij,am->abij');
%     d9 = d9 - permute(d9,[2,1,3,4]) - permute(d9,[1,2,4,3]) + permute(d9,[2,1,4,3]);
%     
%     MM12C = d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9;