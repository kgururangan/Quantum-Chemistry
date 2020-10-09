function [t2a] = update_t2a_ccsdt3(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift)

    % <ijab|H(CCS)|0>
    I2A_vooo = sys.vA_vooo ...
              -0.5*einsum_kg(sys.vA_oooo,t1a,'nmij,an->amij')...
              +einsum_kg(einsum_kg(sys.vA_vovv,t1a,'amef,ei->amif'),t1a,'amif,fj->amij')...
              +einsum_kg(sys.vA_voov,t1a,'amie,ej->amij')...
                   -einsum_kg(sys.vA_voov,t1a,'amje,ei->amij')...
              -0.5*einsum_kg(einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmef,fj->nmej'),t1a,'nmej,an->amej'),t1a,'amej,ei->amij');

    I2A_vvov  = sys.vA_vvov ...
                +0.5*einsum_kg(sys.vA_vvvv,t1a,'abfe,fi->abie')...
                +einsum_kg(einsum_kg(sys.vA_ooov,t1a,'mnie,am->anie'),t1a,'anie,bn->abie');
    
    % H(CCS) intermediates
    h1A_ov = sys.fa_ov...
             +einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me')...
             +einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me'); 
    
    h1A_vv = sys.fa_vv_masked...
             -einsum_kg(h1A_ov,t1a,'me,am->ae')...
             +einsum_kg(sys.vA_vovv,t1a,'amef,fm->ae')...
             +einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae');
     
    h1A_oo = sys.fa_oo_masked...
             +einsum_kg(h1A_ov,t1a,'me,ei->mi')...
             +einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi')...
             +einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi');
     
    h2A_oooo = sys.vA_oooo...
                +einsum_kg(sys.vA_oovo,t1a,'mnej,ei->mnij')...
	     			-einsum_kg(sys.vA_oovo,t1a,'mnei,ej->mnij')...
                +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,ei->mnif'),t1a,'mnif,fj->mnij');
            
    h2A_vvvv = sys.vA_vvvv...
                -einsum_kg(sys.vA_ovvv,t1a,'mbef,am->abef')...
	     			+einsum_kg(sys.vA_ovvv,t1a,'maef,bm->abef')...
                +einsum_kg(einsum_kg(sys.vA_oovv,t1a,'mnef,bn->mbef'),t1a,'mbef,am->abef');
            
    h2A_voov = sys.vA_voov ...
               -einsum_kg(sys.vA_ooov,t1a,'nmie,an->amie')...
               +einsum_kg(sys.vA_vovv,t1a,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vA_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie');
            
    h2B_voov = sys.vB_voov ...
               -einsum_kg(sys.vB_ooov,t1a,'nmie,an->amie')...
               +einsum_kg(sys.vB_vovv,t1a,'amfe,fi->amie')...
               -einsum_kg(einsum_kg(sys.vB_oovv,t1a,'nmfe,fi->nmie'),t1a,'nmie,an->amie');
    
    % VT2 intermediates
    VT2_1A_oo = 0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efin->mi')...
                +einsum_kg(sys.vB_oovv,t2b,'mnef,efin->mi');
            
    VT2_1A_vv = -0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,afmn->ae')...
                -einsum_kg(sys.vB_oovv,t2b,'mnef,afmn->ae');
            
    % note weights of 0.5 and 1 to account on (V T2A^2) and (V T2A T2B),
    % respectively to account for A(ij)A(ab) in contraction term
    VT2_2A_voov = 0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,afin->amie')...
                  +einsum_kg(sys.vB_oovv,t2b,'mnef,afin->amie');  
    
    % weight of 0.5 on (V T2B^2) term to correct for A(ij)A(ab) in
    % contraction term
    VT2_2B_voov = 0.5*einsum_kg(sys.vC_oovv,t2b,'mnef,afin->amie');
    
    VT2_2A_oooo = 0.5*einsum_kg(sys.vA_oovv,t2a,'mnef,efij->mnij');
    
    % contraction intermdiates
    I1A_oo = h1A_oo + VT2_1A_oo;
    I1A_vv = h1A_vv + VT2_1A_vv;
    I2A_voov = h2A_voov + VT2_2A_voov;
    I2B_voov = h2B_voov + VT2_2B_voov;
    I2A_oooo = h2A_oooo + VT2_2A_oooo;
    
    % <ijab| H(CCS) | 0> + <ijab| (H(CCS)T2)_C |0> + <ijab| 0.5*(H(CCS)T2^2)_C |0>
    D1 = -einsum_kg(I2A_vooo,t1a,'amij,bm->abij'); 
    
    D2 = einsum_kg(I2A_vvov,t1a,'abie,ej->abij');
    
    D3 = einsum_kg(I1A_vv,t2a,'ae,ebij->abij');
    
    D4 = -einsum_kg(I1A_oo,t2a,'mi,abmj->abij');
    
    D5 = einsum_kg(I2A_voov,t2a,'amie,ebmj->abij');
    
    D6 = einsum_kg(I2B_voov,t2b,'amie,bejm->abij');
    
    D7 = 0.5*einsum_kg(h2A_vvvv,t2a,'abef,efij->abij');
    
    D8 = 0.5*einsum_kg(I2A_oooo,t2a,'mnij,abmn->abij');
    
    % diagrams that have A(ab)
    D13 = D1 + D3;
    D13 = D13 - permute(D13,[2,1,3,4]);
    
    % diagrams that have A(ij)
    D24 = D2 + D4;
    D24 = D24 - permute(D24,[1,2,4,3]);
        
    % diagrams that have A(ab)A(ij)
    D56 = D5 + D6;
    D56 = D56 - permute(D56,[2,1,3,4]) - permute(D56,[1,2,4,3]) + permute(D56,[2,1,4,3]);
    
    % total CCSD contribution
    X2A = sys.vA_vvoo + D13 + D24 + D56 + D7 + D8;
    
    % CCSD update
    for i = 1:sys.Nocc_alpha
        for j = i+1:sys.Nocc_alpha
            for a = 1:sys.Nvir_alpha
                for b = a+1:sys.Nvir_alpha
                    t2a(a,b,i,j) = X2A(a,b,i,j)/(sys.fa_oo(i,i)+sys.fa_oo(j,j)-sys.fa_vv(a,a)-sys.fa_vv(b,b)-shift);                
                    t2a(b,a,i,j) = -t2a(a,b,i,j);
                    t2a(a,b,j,i) = -t2a(a,b,i,j);
                    t2a(b,a,j,i) = t2a(a,b,i,j);
                end
            end
        end
    end
    
    % CCSDT part    
    PA = sys.PA; PB = sys.PB; HA = sys.HA; HB = sys.HB;
    pA = sys.pA; pB = sys.pB; hA = sys.hA; hB = sys.hB;
    
    h1B_ov = sys.fb_ov + einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me') + einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
    h2B_ooov = sys.vB_ooov + einsum_kg(sys.vB_oovv,t1a,'mnfe,fi->mnie');
    h2A_ooov = sys.vA_ooov + einsum_kg(sys.vA_oovv,t1a,'mnfe,fi->mnie');
    h2A_vovv = sys.vA_vovv - einsum_kg(sys.vA_oovv,t1a,'mnfe,an->amef');
    h2B_vovv = sys.vB_vovv - einsum_kg(sys.vB_oovv,t1a,'nmef,an->amef');
    
    % IJAB
    D1 = einsum_kg(h1A_ov(HA,PA),t3a,'me,abeijm->abij');
    D2 = einsum_kg(h1B_ov(HB,PB),t3b,'me,abeijm->abij');
    D3 = -einsum_kg(h2B_ooov(HA,HB,HA,PB),t3b,'mnif,abfmjn->abij');
    D4 = -0.5*einsum_kg(h2A_ooov(HA,HA,HA,PA),t3a,'mnif,abfmjn->abij');
    D5 = 0.5*einsum_kg(h2A_vovv(PA,HA,PA,PA),t3a,'anef,ebfijn->abij');
    D6 = einsum_kg(h2B_vovv(PA,HB,PA,PB),t3b,'anef,ebfijn->abij');
    D34 = D3 + D4;
    D34 = D34 - permute(D34,[1,2,4,3]);
    D56 = D5 + D6;
    D56 = D56 - permute(D56,[2,1,3,4]);
    X2A_IJAB = D1 + D2 + D34 + D56;
    
    t2a_ABIJ = zeros(sys.Nact_p_alpha,sys.Nact_p_alpha,sys.Nact_h_alpha,sys.Nact_h_alpha);
    for i = 1:sys.Nact_h_alpha
       for j = i+1:sys.Nact_h_alpha
           for a = 1:sys.Nact_p_alpha
               for b = a+1:sys.Nact_p_alpha
                    t2a_ABIJ(a,b,i,j) = X2A_IJAB(a,b,i,j)/(sys.fa_HH(i,i)+sys.fa_HH(j,j)-sys.fa_PP(a,a)-sys.fa_PP(b,b)-shift);                
                    t2a_ABIJ(b,a,i,j) = -t2a_ABIJ(a,b,i,j);
                    t2a_ABIJ(a,b,j,i) = -t2a_ABIJ(a,b,i,j);
                    t2a_ABIJ(b,a,j,i) = t2a_ABIJ(a,b,i,j);
               end
           end
       end
    end
    
    t2a(PA,PA,HA,HA) = t2a(PA,PA,HA,HA) + t2a_ABIJ;
    
    % iJAB
    D1 = -einsum_kg(h2B_ooov(HA,HB,hA,PB),t3b,'mnif,abfmjn->abij');
    D2 = -0.5*einsum_kg(h2A_ooov(HA,HA,hA,PA),t3a,'mnif,abfmjn->abij');
    X2A_iJAB = D1 + D2;
    
    t2a_ABiJ = zeros(sys.Nact_p_alpha,sys.Nact_p_alpha,sys.Nunact_h_alpha,sys.Nact_h_alpha);
    for i = 1:sys.Nunact_h_alpha
        for j = 1:sys.Nact_h_alpha
            for a = 1:sys.Nact_p_alpha
                for b = a+1:sys.Nact_p_alpha
                    t2a_ABiJ(a,b,i,j) = X2A_iJAB(a,b,i,j)/(sys.fa_hh(i,i)+sys.fa_HH(j,j)-sys.fa_PP(a,a)-sys.fa_PP(b,b)-shift); 
                    t2a_ABiJ(b,a,i,j) = -t2a_ABiJ(a,b,i,j);
                end
            end
        end
    end
    
    t2a(PA,PA,hA,HA) = t2a(PA,PA,hA,HA) + t2a_ABiJ;
    
    % IJAb
    D1 = 0.5*einsum_kg(h2A_vovv(pA,HA,PA,PA),t3a,'bnef,aefijn->abij');
    D2 = einsum_kg(h2B_vovv(pA,HB,PA,PB),t3b,'bnef,aefijn->abij');
    X2A_IJAb = D1 + D2;
    
    t2a_AbIJ = zeros(sys.Nact_p_alpha,sys.Nunact_p_alpha,sys.Nact_h_alpha,sys.Nact_h_alpha);
    for i = 1:sys.Nact_h_alpha
        for j = i+1:sys.Nact_h_alpha
            for a = 1:sys.Nact_p_alpha
                for b = 1:sys.Nunact_p_alpha
                    t2a_AbIJ(a,b,i,j) = X2A_IJAb(a,b,i,j)/(sys.fa_HH(i,i)+sys.fa_HH(j,j)-sys.fa_PP(a,a)-sys.fa_pp(b,b)-shift); 
                    t2a_AbIJ(a,b,j,i) = -t2a_AbIJ(a,b,i,j);
                end
            end
        end
    end

    t2a(PA,pA,HA,HA) = t2a(PA,pA,HA,HA) + t2a_AbIJ;



end

