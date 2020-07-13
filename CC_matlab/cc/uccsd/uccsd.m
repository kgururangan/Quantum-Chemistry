function [t1a,t1b,t2a,t2b,t2c,Ecorr] = uccsd(sys,opts)

    diis_size = opts.diis_size;
    maxit = opts.maxit;
    tol = opts.tol;
    
    tic_Start = tic;
    fprintf('\n==================================++Entering UCCSD Routine++=================================\n')

    % Initialize T1 and T2 DIIS containers
    T = zeros(sys.doubles_dim,1);
    T_list = zeros(sys.doubles_dim,diis_size);
    T_resid_list = zeros(sys.doubles_dim,diis_size);

    szt1a = [sys.Nvir_alpha, sys.Nocc_alpha];
    szt1b = [sys.Nvir_beta, sys.Nocc_beta];
    szt2a = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nocc_alpha, sys.Nocc_alpha];
    szt2b = [sys.Nvir_alpha, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_beta];
    szt2c = [sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_beta, sys.Nocc_beta];
    
    % Jacobi/DIIS iterations
    it = 1; flag_conv = 0;
    while it < maxit
        
        tic
        
        % store old T and get current diis dimensions
        T_old = T;
        t1a = reshape(T(sys.posv{1}),szt1a);
        t1b = reshape(T(sys.posv{2}),szt1b);
        t2a = reshape(T(sys.posv{3}),szt2a);
        t2b = reshape(T(sys.posv{4}),szt2b);
        t2c = reshape(T(sys.posv{5}),szt2c);

        % CC correlation energy
        Ecorr = ucc_energy(t1a,t1b,t2a,t2b,t2c,sys);
       
        % update t1 and t2 by Jacobi                        
        [chi1A, chi1B, chi2A, chi2B, chi2C] = build_ucc_intermediates_v2(t1a, t1b, t2a, t2b, t2c, sys);
        t1a = update_t1a(t1a,t1b,t2a,t2b,t2c,chi1A,chi1B,chi2A,chi2B,chi2C,sys);
        t1b = t1a;
        %t2a = update_t2a(t1a,t1b,t2a,t2b,t2c,chi1A,chi1B,chi2A,chi2B,chi2C,sys);
        t2a = update_t2a_jun(t1a,t1b,t2a,t2b,t2c,sys);
        t2b = update_t2b(t1a,t1b,t2a,t2b,t2c,chi1A,chi1B,chi2A,chi2B,chi2C,sys);
        t2c = t2a;
        
        % store vectorized results
        T(sys.posv{1}) = t1a(:); T(sys.posv{2}) = t1b(:);
        T(sys.posv{3}) = t2a(:); T(sys.posv{4}) = t2b(:); T(sys.posv{5}) = t2c(:);

        % build DIIS residual
        T_resid = T - T_old;

        % check for exit condition
        ccsd_resid = sqrt(mean(T_resid.^2));
        if ccsd_resid < tol
            flag_conv = 1;
            break;
        end

        % append trial and residual vectors to lists
        T_list(:,mod(it,diis_size)+1) = T;
        T_resid_list(:,mod(it,diis_size)+1) = T_resid;
         
        % diis extrapolate
        if it >= 2*diis_size
           T = diis_xtrap(T_list,T_resid_list);
        end        

        fprintf('\nIter-%d     Residuum = %4.12f      Ecorr = %4.12f      Elapsed Time = %4.2f s',it,ccsd_resid,Ecorr,toc);
        
        it = it + 1;
%         
    end

    % return final amplitudes and correlation energy
    t1a = reshape(T(sys.posv{1}),szt1a);
    t1b = reshape(T(sys.posv{2}),szt1b);
    t2a = reshape(T(sys.posv{3}),szt2a);
    t2b = reshape(T(sys.posv{4}),szt2b);
    t2c = reshape(T(sys.posv{5}),szt2c);
    Ecorr = ucc_energy(t1a,t1b,t2a,t2b,t2c,sys);
    
    if flag_conv == 1
        fprintf('\nUCCSD successfully converged in %d iterations (%4.2f seconds)\n',it,toc(tic_Start));
        fprintf('Total Energy = %4.12f Ha     Ecorr = %4.12f Ha\n',Ecorr+sys.Escf,Ecorr);
    else
        fprintf('\nUCCSD failed to converged in %d iterations\n',maxit)
    end   
        

end

