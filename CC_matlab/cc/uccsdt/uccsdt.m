function [cc_t,Ecorr] = uccsdt(sys,opts,T_guess)

    diis_size = opts.diis_size;
    maxit = opts.maxit;
    tol = opts.tol;
    shift = opts.shift;

    if nargin < 3 || isempty(T_guess)
        T = zeros(sys.triples_dim,1);
    else
        T = T_guess;
    end
    
    tic_Start = tic;
    fprintf('\n==================================++Entering UCCSDT Routine++=================================\n')

    % Initialize T1, T2, and T3 DIIS containers
    
    T_list = zeros(sys.triples_dim,diis_size);
    T_resid_list = zeros(sys.triples_dim,diis_size);

    szt1a = [sys.Nvir_alpha, sys.Nocc_alpha];
    szt1b = [sys.Nvir_beta, sys.Nocc_beta];
    szt2a = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nocc_alpha, sys.Nocc_alpha];
    szt2b = [sys.Nvir_alpha, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_beta];
    szt2c = [sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_beta, sys.Nocc_beta];
    szt3a = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nvir_alpha, sys.Nocc_alpha, sys.Nocc_alpha, sys.Nocc_alpha];
    szt3b = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_alpha, sys.Nocc_beta];
    szt3c = [sys.Nvir_alpha, sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_beta, sys.Nocc_beta];
    szt3d = [sys.Nvir_beta, sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_beta, sys.Nocc_beta, sys.Nocc_beta];
    
    % Jacobi/DIIS iterations
    it_micro = 0; flag_conv = 0; it_macro = 0;
    while it_micro < maxit
        
        tic
        
        % store old T and get current diis dimensions
        T_old = T;
        
        t1a = reshape(T(sys.posv{1}),szt1a);
        t1b = reshape(T(sys.posv{2}),szt1b);
        t2a = reshape(T(sys.posv{3}),szt2a);
        t2b = reshape(T(sys.posv{4}),szt2b);
        t2c = reshape(T(sys.posv{5}),szt2c);
        t3a = reshape(T(sys.posv{6}),szt3a);
        t3b = reshape(T(sys.posv{7}),szt3b);
        t3c = reshape(T(sys.posv{8}),szt3c);
        t3d = reshape(T(sys.posv{9}),szt3d);

        % CC correlation energy
        Ecorr = ucc_energy(t1a,t1b,t2a,t2b,t2c,sys);
       
        % update t1 and t2 by Jacobi                        
        [HBar_t, VT3_t] = build_ucc_hbar_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);
        
        t1a = update_t1a_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,sys,shift);
        t1b = update_t1b_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,sys,shift);
        t2a = update_t2a_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,sys,shift);
        t2b = update_t2b_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,sys,shift);
        t2c = t2a;
        
        t3a = update_t3a(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,VT3_t,sys,shift);
        t3b = update_t3b(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,VT3_t,sys,shift);
        %t3c = update_t3c(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,VT3_t,sys,shift);
        %t3d = update_t3d(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,HBar_t,VT3_t,sys,shift);
        t3c = permute(t3b,[3,1,2,6,4,5]);
        t3d = t3a;
        
        % store vectorized results
        T(sys.posv{1}) = t1a(:); T(sys.posv{2}) = t1b(:);
        T(sys.posv{3}) = t2a(:); T(sys.posv{4}) = t2b(:); T(sys.posv{5}) = t2c(:);
        T(sys.posv{6}) = t3a(:); T(sys.posv{7}) = t3b(:); T(sys.posv{8}) = t3c(:); T(sys.posv{9}) = t3d(:);

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
        if mod(it_micro,diis_size) == 0 && it_micro > 1
           it_macro = it_macro + 1;
           fprintf('\nDIIS Cycle - %d',it_macro)
           T = diis_xtrap(T_list,T_resid_list);
        end        

        fprintf('\n    Iter-%d     Residuum = %4.12f      Ecorr = %4.12f      Elapsed Time = %4.2f s',it_micro,ccsd_resid,Ecorr,toc);
        
        it_micro = it_micro + 1;
         
    end

    % return final amplitudes and correlation energy
    cc_t.t1a = reshape(T(sys.posv{1}),szt1a);
    cc_t.t1b = reshape(T(sys.posv{2}),szt1b);
    cc_t.t2a = reshape(T(sys.posv{3}),szt2a);
    cc_t.t2b = reshape(T(sys.posv{4}),szt2b);
    cc_t.t2c = reshape(T(sys.posv{5}),szt2c);
    cc_t.t3a = reshape(T(sys.posv{6}),szt3a);
    cc_t.t3b = reshape(T(sys.posv{7}),szt3b);
    cc_t.t3c = reshape(T(sys.posv{8}),szt3c);
    cc_t.t3d = reshape(T(sys.posv{9}),szt3d);
        
    Ecorr = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);
    
    if flag_conv == 1
        fprintf('\nUCCSDT successfully converged in %d iterations (%4.2f seconds)\n',it_micro,toc(tic_Start));
        fprintf('Total Energy = %4.12f Eh     Ecorr = %4.12f Eh\n',Ecorr+sys.Escf,Ecorr);
    else
        fprintf('\nUCCSDT failed to converged in %d iterations\n',maxit)
    end   
        

end

