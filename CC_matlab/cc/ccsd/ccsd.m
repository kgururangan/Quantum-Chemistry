function [t1,t2,Ecorr] = ccsd(sys,opts)

    diis_size = opts.diis_size;
    maxit = opts.maxit;
    tol = opts.tol;
    
    tic_Start = tic;
    
    fprintf('\n==================================++Entering CCSD Routine++=================================\n')
    
    
    % Initialize T1 and T2 DIIS containers
    T = zeros(sys.doubles_dim,1);
    T_list = zeros(sys.doubles_dim,diis_size);
    T_resid_list = zeros(sys.doubles_dim,diis_size);

    % Jacobi/DIIS iterations
    it_micro = 0; flag_conv = 0; it_macro = 0;
    while it_micro < maxit
        
        tic
        
        % store old T and get current diis dimensions
        T_old = T;
        t1 = reshape(T(sys.posv1),sys.Nunocc,sys.Nocc);
        t2 = reshape(T(sys.posv2),sys.Nunocc,sys.Nunocc,sys.Nocc,sys.Nocc);
        
        % CC correlation energy
        Ecorr = cc_energy(t1,t2,sys);
                   
        % update t1 and t2 by Jacobi                        
        t1 = update_t1_ccsd(t1,t2,sys);
        t2 = update_t2_ccsd(t1,t2,sys);
        
        % store vectorized results
        T(sys.posv1) = t1(:); T(sys.posv2) = t2(:);

        % build DIIS residual
        T_resid = T - T_old;
        
        % check for exit condition
        ccsd_resid = sqrt(mean(T_resid.^2));
        if ccsd_resid < tol
            flag_conv = 1;
            break;
        end
        
        % append trial and residual vectors to lists
        T_list(:,mod(it_micro,diis_size)+1) = T;
        T_resid_list(:,mod(it_micro,diis_size)+1) = T_resid;
         
        % diis extrapolate
        if mod(it_micro,diis_size) == 0
           it_macro = it_macro + 1;
           fprintf('\nDIIS Cycle - %d',it_macro)
           T = diis_xtrap(T_list,T_resid_list);
        end        

        fprintf('\n    Iter-%d     Residuum = %4.12f      Ecorr = %4.12f      Elapsed Time = %4.2f s',it_micro,ccsd_resid,Ecorr,toc);
        
        it_micro = it_micro + 1;
        
    end
    
    % return final amplitudes and correlation energy
    t1 = reshape(T(sys.posv1),sys.Nunocc,sys.Nocc);
    t2 = reshape(T(sys.posv2),sys.Nunocc,sys.Nunocc,sys.Nocc,sys.Nocc);
    Ecorr = cc_energy(t1,t2,sys);
    
    if flag_conv == 1
        fprintf('\nCCSD successfully converged in %d iterations (%4.2f seconds)\n',it_micro,toc(tic_Start));
        fprintf('Total Energy = %4.12f Ha     Ecorr = %4.12f Ha\n',Ecorr+sys.Escf,Ecorr);
    else
        fprintf('\nCCSD failed to converged in %d iterations\n',maxit)
    end
            
        

end


       