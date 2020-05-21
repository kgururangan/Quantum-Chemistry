function [L1,L2,lccsd_resid] = lccsd(omega,R,t1,t2,HBar,sys,opts)

% we're using abij ordering for L1 and L2 now!

    fprintf('\n==============++Entering Left-CCSD Routine++=================\n')
    
    FM = sys.FM; VM = sys.VM; occ = sys.occ; unocc = sys.unocc;
    
    tic

    tol = opts.tol;
    diis_size = opts.diis_size;
    maxit = opts.maxit;
    
    Nocc = length(occ); Nunocc = length(unocc);
    Nov = Nocc*Nunocc; Noovv = Nov^2;
    HBar_dim = Nov + Noovv;
    
    LAMBDA_list = zeros(HBar_dim, diis_size); 
    LAMBDA_resid_list = zeros(HBar_dim, diis_size);
    
    if omega == 0.0
        LAMBDA = cat(1,t1(:),t2(:));
        flag_ground = 1;
    else
        LAMBDA = R;
        flag_ground = 0;
    end


    it = 0; flag_conv = 0; 
    while it < maxit
        
        % get current L1 and L2
        L1 = reshape(LAMBDA(1:Nov),Nunocc,Nocc);
        L2 = reshape(LAMBDA(Nov+1:end),Nunocc,Nunocc,Nocc,Nocc);
        
        % build LH diagrammatically with MP denominator separated out
        flag_jacobi = 1;
        [X_ia, X_ijab] = build_L_HBar(L1, L2, HBar, flag_ground, flag_jacobi);
%         
%         X_ia = X_ia - omega*L1;
%         X_ijab = X_ijab - omega*L2;

        % update L1 and L2 by Jacobi
        [L1,L2] = update_l1_l2(X_ia, X_ijab, omega, HBar, occ, unocc);
        LAMBDA = cat(1,L1(:),L2(:));
        
        % buid LH - omega*L residual measure (use full LH)
        flag_jacobi = 0;
        [LH1,LH2] = build_L_HBar(L1, L2, HBar, flag_ground, flag_jacobi);
        LH = cat(1,LH1(:),LH2(:));
        LAMBDA_resid = LH - omega*LAMBDA;
        
        % get lcc energy - returns Ecorr + omega
        E_lcc = lcc_energy(LH,t1,t2,VM,FM,occ,unocc);
        
        % check exit condition
        lccsd_resid = norm(LAMBDA_resid);
        if lccsd_resid < tol
            flag_conv = 1;
            break;
        end
        
        % append trial and residual vectors to diis lists
        LAMBDA_list(:,mod(it,diis_size)+1) = LAMBDA;
        LAMBDA_resid_list(:,mod(it,diis_size)+1) = LAMBDA_resid;
 
        % diis extrapolate
        if it >= diis_size
           LAMBDA = diis_xtrap(LAMBDA_list,LAMBDA_resid_list);
        end
        
        % biorthonormalize to right eigenvectors
        % note |R| = 1, so this necesarily enforces |L| = 1
        if flag_ground ~= 1
            LAMBDA = LAMBDA./(LAMBDA'*R);
        end

        % Print status
        fprintf('\nIter-%d        Energy = %4.8f      Residuum = %4.8f',...
                    it, E_lcc, lccsd_resid);
        
        it = it + 1;
        
    end
    
    if flag_conv == 1
        fprintf('\nLeft-CCSD successfully converged in %d iterations (%4.2f seconds)\n',it,toc);
    else
        fprintf('\nLeft-CCSD failed to converged in %d iterations\n',maxit)
    end
         

end

    
  