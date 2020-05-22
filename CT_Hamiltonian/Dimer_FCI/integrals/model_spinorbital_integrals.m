function [Zmat, Vmat] = model_spinorbital_integrals(par, n_sp_orb)

    J_S_S = par.singlet_singlet;
    J_T_T = par.triplet_triplet;
    J_CT_T = par.ct_triplet;
    J_CT_S = par.ct_singlet;
    J_CT_CT = par.ct_ct;
    J_S_T = par.singlet_triplet;
    
    for p = 1:n_sp_orb
        for q = 1:n_sp_orb
            
            if mod(p,2) == mod(q,2)
                Zmat(p,q) = 
            
            
            
            for r = 1:n_sp_orb
                for s = 1:n_sp_orb
                    
                end
            end
        end
    end
    

end

