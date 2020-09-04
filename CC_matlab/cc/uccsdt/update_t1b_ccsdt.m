function [t1b] = update_t1b_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,sys,shift)

    chi1B_vv = sys.fb_vv_masked + einsum_kg(sys.vC_vovv,t1b,'anef,fn->ae') + einsum_kg(sys.vB_ovvv,t1a,'nafe,fn->ae');
    chi1B_oo = sys.fb_oo_masked + einsum_kg(sys.vC_ooov,t1b,'mnif,fn->mi') + einsum_kg(sys.vB_oovo,t1b,'nmfi,fn->mi');
    h1A_ov = sys.fa_ov + einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me') + einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');
    h1B_ov = sys.fb_ov + einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me') + einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
    h1B_oo = chi1B_oo + einsum_kg(h1B_ov,t1b,'me,ei->mi');
       
    M11 = sys.fb_vo - einsum_kg(h1B_oo,t1b,'mi,am->ai') + einsum_kg(chi1B_vv,t1b,'ae,ei->ai')...
          +einsum_kg(sys.vC_voov,t1b,'anif,fn->ai') + einsum_kg(sys.vB_ovvo,t1a,'nafi,fn->ai');
      
    h2C_ooov = sys.vC_ooov + einsum_kg(sys.vC_oovv,t1b,'mnfe,fi->mnie'); 
    h2B_oovo = sys.vB_oovo + einsum_kg(sys.vB_oovv,t1b,'nmef,fi->nmei');
    h2C_vovv = sys.vC_vovv - einsum_kg(sys.vC_oovv,t1b,'mnfe,an->amef');
    h2B_ovvv = sys.vB_ovvv - einsum_kg(sys.vB_oovv,t1b,'mnfe,an->mafe');
      
    CCS_T2 = +einsum_kg(h1A_ov,t2b,'me,eami->ai') + einsum_kg(h1B_ov,t2c,'me,aeim->ai')...
             -0.5*einsum_kg(h2C_ooov,t2c,'mnif,afmn->ai') - einsum_kg(h2B_oovo,t2b,'nmfi,fanm->ai')...
             +0.5*einsum_kg(h2C_vovv,t2c,'anef,efin->ai') + einsum_kg(h2B_ovvv,t2b,'nafe,feni->ai');
       
    X1B = M11 + CCS_T2; 
    
    % CCSDT part
    X1B = X1B ...
          + 0.25*einsum_kg(sys.vC_oovv,t3d,'mnef,aefimn->ai')...
          + 0.25*einsum_kg(sys.vA_oovv,t3b,'mnef,efamni->ai')...
          + einsum_kg(sys.vB_oovv,t3c,'mnef,efamni->ai');
    
    for a = 1:sys.Nvir_beta
        for i = 1:sys.Nocc_beta
            denom = sys.fb_oo(i,i) - sys.fb_vv(a,a) - shift;
            t1b(a,i) = X1B(a,i) / denom;
        end
    end


end