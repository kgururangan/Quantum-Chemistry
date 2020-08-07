function [cc_t,lccsd_resid] = luccsd(cc_t,HBar_t,sys,opts)

    fprintf('\n==================================++Entering Left-UCCSD Routine++=============================\n')
    
    tic_Start = tic;

    tol = opts.tol;
    diis_size = opts.diis_size;
    maxit = opts.maxit;
    shift = opts.shift;
    
    Ecorr = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);
    
    % Initialize T1 and T2 DIIS containers
    LAMBDA = cat(1,cc_t.t1a(:),cc_t.t1b(:),cc_t.t2a(:),cc_t.t2b(:),cc_t.t2c(:));
    LAMBDA_list = zeros(sys.doubles_dim,diis_size);
    LAMBDA_resid_list = zeros(sys.doubles_dim,diis_size);

    szl1a = [sys.Nvir_alpha, sys.Nocc_alpha];
    szl1b = [sys.Nvir_beta, sys.Nocc_beta];
    szl2a = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nocc_alpha, sys.Nocc_alpha];
    szl2b = [sys.Nvir_alpha, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_beta];
    szl2c = [sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_beta, sys.Nocc_beta];
    
    omega = 0.0;
    flag_ground = true;

    it_micro = 0; it_macro = 0; flag_conv = 0; 
    while it_micro < maxit

        tic
        
        % get current L1 and L2
        l1a = reshape(LAMBDA(sys.posv{1}),szl1a);
        l1b = reshape(LAMBDA(sys.posv{2}),szl1b);
        l2a = reshape(LAMBDA(sys.posv{3}),szl2a);
        l2b = reshape(LAMBDA(sys.posv{4}),szl2b);
        l2c = reshape(LAMBDA(sys.posv{5}),szl2c);

        % update L1 and L2 by Jacobi
        flag_jacobi = true;
        X1A = build_LH_1A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        X1B = build_LH_1B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        X2A = build_LH_2A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        X2B = build_LH_2B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        X2C = build_LH_2C(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        
        [l1a, l1b, l2a, l2b, l2c] = update_L(X1A,X1B,X2A,X2B,X2C,HBar_t,sys,omega,shift);
       
        LAMBDA = cat(1,l1a(:),l1b(:),l2a(:),l2b(:),l2c(:));
        
        % buid LH - omega*L residual measure (use full LH)
        [LH] = build_ucc_LHBar(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_ground);
        LAMBDA_resid = LH - omega*LAMBDA;

        % get lcc energy - returns Ecorr + omega
        E_lcc = sqrt(sum(LH.^2))./sqrt(sum(LAMBDA.^2)) + Ecorr;
        
        % check exit condition
        lccsd_resid = sqrt(sum(LAMBDA_resid.^2));
        if lccsd_resid < tol
            flag_conv = 1;
            break;
        end
        
        % append trial and residual vectors to diis lists
        LAMBDA_list(:,mod(it_micro,diis_size)+1) = LAMBDA;
        LAMBDA_resid_list(:,mod(it_micro,diis_size)+1) = LAMBDA_resid;
 
        % diis extrapolate
        if mod(it_micro,diis_size) == 0 && it_micro > 1
           it_macro = it_macro + 1;
           fprintf('\nDIIS Cycle - %d',it_macro)
           LAMBDA = diis_xtrap(LAMBDA_list,LAMBDA_resid_list);
        end

        % Print status
        fprintf('\nIter-%d     Residuum = %4.12f      Ecorr = %4.12f      Elapsed Time = %4.2f s',it_micro,lccsd_resid,E_lcc,toc);
        
        it_micro = it_micro + 1;
        
    end
    
    cc_t.l1a{1} = l1a;
    cc_t.l1b{1} = l1b;
    cc_t.l2a{1} = l2a;
    cc_t.l2b{1} = l2b;
    cc_t.l2c{1} = l2c;
    
    if flag_conv == 1
        fprintf('\nLeft-UCCSD successfully converged in %d iterations (%4.2f seconds)\n',it_micro,toc(tic_Start));
    else
        fprintf('\nLeft-UCCSD failed to converged in %d iterations\n',maxit)
    end
         

end