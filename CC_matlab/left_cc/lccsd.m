function [L1,L2,lccsd_resid] = lccsd(t1,t2,HBar,sys,opts)

% we're using abij ordering for L1 and L2 now!

    fprintf('\n==================================++Entering Left-CCSD Routine++=============================\n')
    
    %FM = sys.FM; VM = sys.VM; occ = sys.occ; unocc = sys.unocc;
    
    tic_Start = tic;

    tol = opts.tol;
    diis_size = opts.diis_size;
    maxit = opts.maxit;
    shift = opts.shift;
    Nocc = sys.Nocc; Nunocc = sys.Nunocc;
    doubles_dim = sys.doubles_dim;
    posv1 = sys.posv1;
    posv2 = sys.posv2;
    
    size_L1 = [Nunocc, Nocc];
    size_L2 = [Nunocc, Nunocc, Nocc, Nocc];
    
    LAMBDA_list = zeros(doubles_dim, diis_size); 
    LAMBDA_resid_list = zeros(doubles_dim, diis_size);
    
    LAMBDA = cat(1,t1(:),t2(:));
    flag_ground = 1;
    omega = 0.0;


    it_micro = 0; it_macro = 0; flag_conv = 0; 
    while it_micro < maxit

        tic
        
        % get current L1 and L2
        L1 = reshape(LAMBDA(posv1),size_L1);
        L2 = reshape(LAMBDA(posv2),size_L2);
        
        % build LH diagrammatically with MP denominator separated out
        flag_jacobi = 1;
        [X_ia, X_ijab] = build_L_HBar(L1, L2, HBar, flag_ground, flag_jacobi);

        % update L1 and L2 by Jacobi
        [L1,L2] = update_l1_l2(X_ia, X_ijab, omega, HBar, shift);
        LAMBDA = cat(1,L1(:),L2(:));
        
        % buid LH - omega*L residual measure (use full LH)
        flag_jacobi = 0;
        [LH1,LH2] = build_L_HBar(L1, L2, HBar, flag_ground, flag_jacobi);
        LH = cat(1,LH1(:),LH2(:));
        LAMBDA_resid = LH;
        
        % get lcc energy - returns Ecorr + omega
        E_lcc = lcc_energy(LAMBDA,LH,t1,t2,sys);
        
        % check exit condition
        lccsd_resid = norm(LAMBDA_resid);
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
    
    if flag_conv == 1
        fprintf('\nLeft-CCSD successfully converged in %d iterations (%4.2f seconds)\n',it_micro,toc(tic_Start));
    else
        fprintf('\nLeft-CCSD failed to converged in %d iterations\n',maxit)
    end
         

end

    
  