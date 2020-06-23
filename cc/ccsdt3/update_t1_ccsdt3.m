function [t1] = update_t1_ccsdt3(t1,t2,t3,sys)

    chi_me = sys.fov + einsum_kg(sys.Voovv,t1,'mnef,fn->me');
    
    chi_ae = sys.fvv_masked + ...
             einsum_kg(sys.Vvovv,t1,'anef,fn->ae') - ...
             0.5*einsum_kg(sys.Voovv,t2,'mnef,afmn->ae') - ...
             einsum_kg(chi_me,t1,'me,am->ae');
         
    chi_mi = sys.foo_masked + ...
             einsum_kg(sys.Vooov,t1,'mnif,fn->mi') + ...
             0.5*einsum_kg(sys.Voovv,t2,'mnef,efin->mi');

    % build T1
    
    % < phi_{ia} | H_{CCSD} | 0 >
    TEMP1 = einsum_kg(chi_me,t2,'me,aeim->ai') - ...
            einsum_kg(sys.Vovov,t1,'maie,em->ai') - ...
            einsum_kg(chi_mi,t1,'mi,am->ai') + ...
            einsum_kg(chi_ae,t1,'ae,ei->ai') - ...
            0.5*einsum_kg(sys.Vooov,t2,'mnif,afmn->ai') + ...
            0.5*einsum_kg(sys.Vvovv,t2,'anef,efin->ai');
        
    % < phi_{ia} | (H_N T3)_C | 0 >
    TEMP1_3 = 0.25*einsum_kg(sys.VHHPP,t3,'mnef,aefimn->ai');
    
    X_ai = sys.fvo + TEMP1; 
    
    for a = 1:sys.Nunocc
        for i = 1:sys.Nocc
            if a <= sys.Nact_p && i - sys.Ncore > 0
                t1(a,i) = (X_ai(a,i)+TEMP1_3(a,i-sys.Ncore))/(sys.foo(i,i)-sys.fvv(a,a));
            else
                t1(a,i) = X_ai(a,i)/(sys.foo(i,i)-sys.fvv(a,a));
            end
        end
    end
    
end
    