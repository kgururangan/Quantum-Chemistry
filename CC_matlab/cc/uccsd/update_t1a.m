function [t1a] = update_t1a(t1a,t1b,t2a,t2b,t2c,sys,shift)

    chi1A_vv = sys.fa_vv + einsum_kg(sys.vA_vovv,t1a,'anef,fn->ae') + einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae');
    chi1A_oo = sys.fa_oo + einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi') + einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi');
    h1A_ov = sys.fa_ov + einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me') + einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');
    h1B_ov = sys.fb_ov + einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me') + einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
    h1A_oo = chi1A_oo + einsum_kg(h1A_ov,t1a,'me,ei->mi');
       
    M11 = sys.fa_vo - einsum_kg(h1A_oo,t1a,'mi,am->ai') + einsum_kg(chi1A_vv,t1a,'ae,ei->ai')...
          +einsum_kg(sys.vA_voov,t1a,'anif,fn->ai') + einsum_kg(sys.vB_voov,t1b,'anif,fn->ai');
      
    h2A_ooov = sys.vA_ooov + einsum_kg(sys.vA_oovv,t1a,'mnfe,fi->mnie'); 
    h2B_ooov = sys.vB_ooov + einsum_kg(sys.vB_oovv,t1a,'mnfe,fi->mnie');
    h2A_vovv = sys.vA_vovv - einsum_kg(sys.vA_oovv,t1a,'mnfe,an->amef');
    h2B_vovv = sys.vB_vovv - einsum_kg(sys.vB_oovv,t1a,'nmef,an->amef');
      
    CCS_T2 = +einsum_kg(h1A_ov,t2a,'me,aeim->ai') + einsum_kg(h1B_ov,t2b,'me,aeim->ai')...
             -0.5*einsum_kg(h2A_ooov,t2a,'mnif,afmn->ai') - einsum_kg(h2B_ooov,t2b,'mnif,afmn->ai')...
             +0.5*einsum_kg(h2A_vovv,t2a,'anef,efin->ai') + einsum_kg(h2B_vovv,t2b,'anef,efin->ai');
       
    X1A = M11 + CCS_T2; 
    
    for a = 1:sys.Nvir_alpha
        for i = 1:sys.Nocc_alpha
            denom = sys.fa_oo(i,i) - sys.fa_vv(a,a) - shift;
            t1a(a,i) = t1a(a,i) + X1A(a,i) / denom;
        end
    end


end

