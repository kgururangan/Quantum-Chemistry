function [HBar_t] = build_hbar_ccs(t1a,t1b,sys)

    % h(ov)
    H1A.ov = sys.fa_ov...
             +einsum_kg(sys.vA_oovv,t1a,'mnef,fn->me')...
             +einsum_kg(sys.vB_oovv,t1b,'mnef,fn->me');         
    H1B.ov = sys.fb_ov...
             +einsum_kg(sys.vB_oovv,t1a,'nmfe,fn->me')...
             +einsum_kg(sys.vC_oovv,t1b,'mnef,fn->me');
         
    % h(oo)
    H1A.oo = sys.fa_oo...
             +einsum_kg(H1A.ov,t1a,'me,ei->mi')...
             +einsum_kg(sys.vA_ooov,t1a,'mnif,fn->mi')...
             +einsum_kg(sys.vB_ooov,t1b,'mnif,fn->mi');
    H1B.oo = sys.fb_oo...
             +einsum_kg(H1B.ov,t1b,'me,ei->mi')...
             +einsum_kg(sys.vB_oovo,t1a,'nmfi,fn->mi')...
             +einsum_kg(sys.vC_ooov,t1b,'mnif,fn->mi');
         
    % h(vv)
    H1A.vv = sys.fa_vv...
             -einsum_kg(H1A.ov,t1a,'me,am->ae')...
             +einsum_kg(sys.vA_vovv,t1a,'amef,fm->ae')...
             +einsum_kg(sys.vB_vovv,t1b,'anef,fn->ae');
    H1B.vv = sys.fb_vv...
             -einsum_kg(H1B.ov,t1b,'me,am->ae')...
             +einsum_kg(sys.vB_ovvv,t1a,'nafe,fn->ae')...
             +einsum_kg(sys.vC_vovv,t1b,'anef,fn->ae');
         
    % h(vovv)
    H2A.vovv = sys.vA_vovv - einsum_kg(sys.vA_oovv,t1a,'mnfe,an->amef');
    H2B.vovv = sys.vB_vovv - einsum_kg(sys.vB_oovv,t1a,'nmef,an->amef');
    H2B.ovvv = sys.vB_ovvv - einsum_kg(sys.vB_oovv,t1b,'mnfe,an->mafe');
    H2C.vovv = sys.vC_vovv - einsum_kg(sys.vC_oovv,t1b,'mnfe,an->amef');
    
    % h(ooov)
    H2A.ooov = sys.vA_ooov + einsum_kg(sys.vA_oovv,t1a,'mnfe,fi->mnie');
    H2B.ooov = sys.vB_ooov + einsum_kg(sys.vB_oovv,t1a,'mnfe,fi->mnie');
    H2B.oovo = sys.vB_oovo + einsum_kg(sys.vB_oovv,t1b,'nmef,fi->nmei');
    H2C.ooov = sys.vC_ooov + einsum_kg(sys.vC_oovv,t1b,'mnfe,fi->mnie');
    
    % h(oooo)
    H2A.oooo = sys.vA_oooo + einsum_kg(H2A.ooov,t1a,'nmje,ei->mnij') - einsum_kg(sys.vA_oovo,t1a,'mnei,ej->mnij'); 
    H2
    
    
    % store HBar into HBar_t struct
    HBar_t.H1A = H1A;
    HBar_t.H1B = H1B;
    HBar_t.H2A = H2A;
    HBar_t.H2B = H2B;
    HBar_t.H2C = H2C;

end

