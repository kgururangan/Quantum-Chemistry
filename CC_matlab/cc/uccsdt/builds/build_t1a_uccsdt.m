function [X1a_ai] = build_t1a_uccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys)

    [chi1A, chi1B, chi2A, chi2B, chi2C] = build_ucc_intermediates_v2(t1a, t1b, t2a, t2b, t2c, sys);
    [HBar_t, VT3_t] = build_ucc_hbar_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);


    % CCSD part
    X1a_ai = sys.fa_vo +... % 1
             +einsum_kg(chi1A.ae_bar,t1a,'ae,ei->ai')... % 8
             -einsum_kg(chi1A.mi_bar,t1a,'mi,am->ai')... % 8 (effectively 5... 3 in common between chi.ae and chi.mi)
             +einsum_kg(chi1A.me,t2a,'me,aeim->ai') ... % 3
             +einsum_kg(chi1B.me,t2b,'me,aeim->ai') ... % 3
             +0.5*einsum_kg(sys.vA_vovv,t2a,'amef,efim->ai')... % 1
             -0.5*einsum_kg(sys.vA_ooov,t2a,'mnie,aemn->ai')... % 1
             +einsum_kg(sys.vB_vovv,t2b,'amef,efim->ai')... % 1
             -einsum_kg(sys.vB_ooov,t2b,'mnie,aemn->ai')... % 1
             +einsum_kg(sys.vA_voov,t1a,'amie,em->ai')... % 1
             +einsum_kg(sys.vB_voov,t1b,'amie,em->ai'); % 1

    % CCSDT part
    X1a_ai = X1a_ai...
             +0.25*einsum_kg(sys.vA_oovv,t3a,'mnef,aefimn->ai')...
             +einsum_kg(sys.vB_oovv,t3b,'mnef,aefimn->ai')...
             +0.25*einsum_kg(sys.vC_oovv,t3c,'mnef,aefimn->ai');
         
end