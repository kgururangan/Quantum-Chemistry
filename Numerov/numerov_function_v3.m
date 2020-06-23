function [ENERGY, PSI, energy_hist, delta_hist, flag_exits] = numerov_function_v3(x,U,mass,Elim,NNc,maxit,tol,maxreit,flag_verbose)


    % Initialize containers to store parameters for each node solution
    ENERGY = []; ENERGY_HIST = {};
    PSI = []; DELTA_HIST = {};
    ct_nodes = 1; reit = 0;
    
    for num_nodes = NNc+1 % accounting for extra node at x = 0 in the node count

        % initialize energy bounds for given node count
        if isempty(Elim)
            EL = 0; EH = U(end);
        else
            EL = Elim(ct_nodes,1); EH = Elim(ct_nodes,2);
        end
        
        [energy, psi, energy_hist, delta_hist, flag_exit] =...
            numerov_function_v2(x,U,mass,[EL,EH],num_nodes,maxit,tol,flag_verbose);
        
        while flag_exit ~= 0 && reit < maxreit
            fprintf('Restarting Numerov with new energy bounds...\n')
            EL = max(energy-0.25*energy,0); EH = min(energy+0.25*energy,U(end));
            EL
            EH
            [energy, psi, energy_hist, delta_hist, flag_exit] =...
                  numerov_function_v2(x,U,mass,[EL,EH],num_nodes,maxit,tol,flag_verbose);
            reit = reit + 1;
        end

        % Record at given node count
        ENERGY(ct_nodes) = energy;
        PSI(ct_nodes,:) = psi;
        ENERGY_HIST{ct_nodes} = energy_hist;
        DELTA_HIST{ct_nodes} = delta_hist;
        flag_exits(ct_nodes) = flag_exit;
        
        ct_nodes = ct_nodes + 1;
    end
    
end
            
