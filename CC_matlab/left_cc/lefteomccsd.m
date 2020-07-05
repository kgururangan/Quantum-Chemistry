function [LAMBDA,E_lcc,lccsd_resid] = lefteomccsd(omega,Rvec,HBar,t1,t2,sys,opts)
% we're using abij ordering for L1 and L2 now!

    fprintf('\n==================================++Entering Left-EOM-CCSD Routine++===========================\n')
    
    tic

    tol = opts.tol;
    diis_size = opts.diis_size;
    maxit = opts.maxit;
    nroot = length(omega);

    Nunocc = sys.Nunocc; Nocc = sys.Nocc;

    LAMBDA_list = cell(1,nroot);
    LAMBDA_resid_list = cell(1,nroot);
    lccsd_resid = zeros(1,nroot);
    E_lcc = zeros(1,nroot);
    
    for j = 1:nroot
        LAMBDA_list{j} = zeros(sys.doubles_dim, diis_size);
        LAMBDA_resid_list{j} = zeros(sys.doubles_dim, diis_size); 
    end

    LAMBDA = Rvec(:,1:nroot);

    idx_not_converged = 1:nroot;

    it = 0; flag_conv = 0; flag_ground = 0;
    while it < maxit && flag_conv == 0
        
        fprintf('\nIter-%d:\n', it)
        fprintf('-------------------------------------------------------\n')
        
        for j = idx_not_converged
            
            % get current L1 and L2
            L1 = reshape(LAMBDA(sys.posv1,j),Nunocc,Nocc);
            L2 = reshape(LAMBDA(sys.posv2,j),Nunocc,Nunocc,Nocc,Nocc);

            % build LH diagrammatically with MP denominator separated out
            flag_jacobi = 1;
            [X_ia, X_ijab] = build_L_HBar(L1, L2, HBar, flag_ground, flag_jacobi);

            % update L1 and L2 by Jacobi
            [L1,L2] = update_l1_l2(X_ia, X_ijab, omega(j), HBar);
            LAMBDA(:,j) = cat(1,L1(:),L2(:));

            % buid LH - omega*L residual measure (use full LH)
            flag_jacobi = 0;
            [LH1,LH2] = build_L_HBar(L1, L2, HBar, flag_ground, flag_jacobi);
            LH = cat(1,LH1(:),LH2(:));
            LAMBDA_resid = LH - omega(j)*LAMBDA(:,j);

            % get lcc energy - returns Ecorr + omega
            E_lcc(j) = lcc_energy(LAMBDA(:,j),LH,t1,t2,sys);

            % check exit condition
            lccsd_resid(j) = norm(LAMBDA_resid);

            % append trial and residual vectors to diis lists
            LAMBDA_list{j}(:,mod(it,diis_size)+1) = LAMBDA(:,j);
            LAMBDA_resid_list{j}(:,mod(it,diis_size)+1) = LAMBDA_resid;

            % diis extrapolate
            if it >= diis_size
               LAMBDA(:,j) = diis_xtrap(LAMBDA_list{j},LAMBDA_resid_list{j});
            end

        end
        
        % biorthonormalize to right eigenvectors using full matrices
        LAMBDA = LAMBDA/(LAMBDA'*Rvec(:,1:nroot));
        
        % print status
        for j = 1:nroot
            fprintf('   Root - %d     e = %4.10f     |r| = %4.10f\n',j,E_lcc(j),lccsd_resid(j));
        end

        idx_not_converged = find(lccsd_resid > tol);
        
        if all(lccsd_resid < tol)
            flag_conv = 1;
        end
        
        it = it + 1;
        
    end
    
    % build final residuals
    for j = 1:nroot
        [L1, L2] = build_L_HBar(reshape(LAMBDA(sys.posv1,j),Nunocc,Nocc),...
                                reshape(LAMBDA(sys.posv2,j),Nunocc,Nunocc,Nocc,Nocc),...
                                HBar, 0, 0);
                    
        LH = cat(1,L1(:),L2(:));
    
        lccsd_resid(j) = norm(LH - omega(j)*LAMBDA(:,j));
        E_lcc(j) = real(lcc_energy(LAMBDA(:,j),LH,t1,t2,sys));
    end
    
    if flag_conv == 1
        fprintf('\nLeft-EOMCCSD successfully converged in %d iterations (%4.2f seconds)\n',it,toc);
    else
        fprintf('\nLeft-EOMCCSD failed to converged in %d iterations\n',maxit)
    end
         
    biorthmat = transpose(conj(LAMBDA))*Rvec(:,1:nroot);
    fprintf('|LR - 1| = %4.12f\n',norm(biorthmat-eye(nroot)))
    
end


