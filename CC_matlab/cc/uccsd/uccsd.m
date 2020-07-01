function [t1a,t1b,t2a,t2b,t2c,Ecorr] = uccsd(sys,opts)

    diis_size = opts.diis_size;
    maxit = opts.maxit;
    tol = opts.tol;
    
    tic_Start = tic;
    fprintf('\n==================================++Entering UCCSD Routine++=================================\n')
    
    t1a = zeros(sys.Nvir_alpha,sys.Nocc_alpha);
    t1b = zeros(sys.Nvir_beta,sys.Nocc_beta);
    t2a = zeros(sys.Nvir_alpha,sys.Nvir_alpha,sys.Nocc_alpha,sys.Nocc_alpha);
    t2b = zeros(sys.Nvir_alpha,sys.Nvir_beta,sys.Nocc_alpha,sys.Nocc_beta);
    t2c = zeros(sys.Nvir_beta,sys.Nvir_beta,sys.Nocc_beta,sys.Nocc_beta);

    
    % Initialize T1 and T2 DIIS containers
    T1A = zeros(numel(t1a),1);
    T1B = zeros(numel(t1b),1);
    T2A = zeros(numel(t2a),1);
    T2B = zeros(numel(t2b),1);
    T2C = zeros(numel(t2c),1);
    
    T1A_list = zeros(numel(t1a),1);
    T1B_list = zeros(numel(t1b),1);
    T2A_list = zeros(numel(t2a),1);
    T2B_list = zeros(numel(t2b),1);
    T2C_list = zeros(numel(t2c),1);
    
    T1A_resid_list = zeros(numel(t1a),1);
    T1B_resid_list = zeros(numel(t1b),1);
    T2A_resid_list = zeros(numel(t2a),1);
    T2B_resid_list = zeros(numel(t2b),1);
    T2C_resid_list = zeros(numel(t2c),1);
    
    % Jacobi/DIIS iterations
    it = 1; flag_conv = 0;
    while it < maxit
        
        tic
        
        % store old T and get current diis dimensions
        T1A_old = T1A;
        T1B_old = T1B;
        T2A_old = T2A;
        T2B_old = T2B;
        T2C_old = T2C;

        curr_size = size(T1A_list,2);
        
        % build t1 and t2 diagrammatically
        [chi1A, chi1B, chi2A, chi2B, chi2C] = build_ucc_intermediates_v2(t1a, t1b, t2a, t2b, t2c, sys);
        t1a = update_t1a(t1a,t1b,t2a,t2b,t2c,chi1A,chi1B,chi2A,chi2B,chi2C,sys);
        t1b = t1a;
        t2a = update_t2a(t1a,t1b,t2a,t2b,t2c,chi1A,chi1B,chi2A,chi2B,chi2C,sys);
        t2b = update_t2b(t1a,t1b,t2a,t2b,t2c,chi1A,chi1B,chi2A,chi2B,chi2C,sys);
        t2c = t2a;

        
        % store vectorized results
        T1A = t1a(:); T1B = t1b(:); T2A = t2a(:); T2B = t2b(:); T2C = t2c(:);

        % build DIIS residual
        T1A_resid = T1A - T1A_old;
        T1B_resid = T1B - T1B_old;
        T2A_resid = T2A - T2A_old;
        T2B_resid = T2B - T2B_old;
        T2C_resid = T2C - T2C_old;
        
        % check for exit condition
        ccsd_resid = sqrt( mean(T1A_resid.^2) +...
                           mean(T1B_resid.^2) +...
                           mean(T2A_resid.^2) +...
                           mean(T2B_resid.^2) +...
                           mean(T2C_resid.^2) );
                       
        if ccsd_resid < tol
            flag_conv = 1;
            break;
        end
        
        % append trial and residual vectors to lists
        T1A_list(:,curr_size+1) = T1A;
        T1B_list(:,curr_size+1) = T1B;
        T2A_list(:,curr_size+1) = T2A;
        T2B_list(:,curr_size+1) = T2B;
        T2C_list(:,curr_size+1) = T2C;
        
        T1A_resid_list(:,curr_size+1) = T1A_resid;
        T1B_resid_list(:,curr_size+1) = T1B_resid;
        T2A_resid_list(:,curr_size+1) = T2A_resid;
        T2B_resid_list(:,curr_size+1) = T2B_resid;
        T2C_resid_list(:,curr_size+1) = T2C_resid;
        
        % enforce max diis_size
        if size(T1A_list,2) > diis_size
            
            T1A_list(:,1) = [];
            T1B_list(:,1) = [];
            T2A_list(:,1) = [];
            T2B_list(:,1) = [];
            T2C_list(:,1) = [];
            
            T1A_resid_list(:,1) = [];
            T1B_resid_list(:,1) = [];
            T2A_resid_list(:,1) = [];
            T2B_resid_list(:,1) = [];
            T2C_resid_list(:,1) = [];
            
        end
        
        % diis extrapolate
        if it >= 2
           T1A = diis_xtrap(T1A_list,T1A_resid_list);
           T1B = diis_xtrap(T1B_list,T1B_resid_list);
           T2A = diis_xtrap(T2A_list,T2A_resid_list);
           T2B = diis_xtrap(T2B_list,T2B_resid_list);
           T2C = diis_xtrap(T2C_list,T2C_resid_list);
        end
        
        t1a = reshape(T1A,size(t1a));
        t1b = reshape(T1B,size(t1b));
        t2a = reshape(T2A,size(t2a));
        t2b = reshape(T2B,size(t2b));
        t2c = reshape(T2C,size(t2c));

        % CC correlation energy
        Ecorr = ucc_energy(t1a,t1b,t2a,t2b,t2c,sys);
        
        fprintf('\nIter-%d     Residuum = %4.12f      Ecorr = %4.12f      Elapsed Time = %4.2f s',it,ccsd_resid,Ecorr,toc);
        
        it = it + 1;
        
    end
    
    t1a = reshape(T1A,size(t1a));
    t1b = reshape(T1B,size(t1b));
    t2a = reshape(T2A,size(t2a));
    t2b = reshape(T2B,size(t2b));
    t2c = reshape(T2C,size(t2c));
    Ecorr = ucc_energy(t1a,t1b,t2a,t2b,t2c,sys);
    
    
    if flag_conv == 1
        fprintf('\nUCCSD successfully converged in %d iterations (%4.2f seconds)\n',it,toc(tic_Start));
        fprintf('Final Correlation Energy = %4.12f Ha\n',Ecorr);
    else
        fprintf('\nUCCSD failed to converged in %d iterations\n',maxit)
    end
            
        

end

