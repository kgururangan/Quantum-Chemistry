function [t1a] = update_t1a(t1a, t1b, t2a, t2b, t2c, chi1A, chi1B, chi2A, chi2B, chi2C, sys, shift)

    X1a_ia = sys.fa_vo +... % 1
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

         
    omega = 1;          
    for a = 1:sys.Nvir_alpha
        for i = 1:sys.Nocc_alpha
            temp = X1a_ia(a,i)/(sys.fa_oo(i,i) - sys.fa_vv(a,a) - shift);
            t1a(a,i) = (1-omega)*t1a(a,i) + omega*temp;
        end
    end
             
         
end