function [LAMBDA,lccsd_resid,cc_t,E_lcc] = leftcc_solver2(omega,Rvec,HBar_t,cc_t,sys,opts)
% we're using abij ordering for L1 and L2 now!

    fprintf('\n==================================++Entering Left-EOM-UCCSD Routine++===========================\n')
    
    tic_Start = tic;

    tol = opts.tol;
    diis_size = opts.diis_size;
    maxit = opts.maxit;
    shift = opts.shift;
    nroot = opts.nroot;

    Ecorr = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);

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

%     idx_root = @(j) [(j-1)*sys.doubles_dim+1:j*sys.doubles_dim];

    it_micro = 0; it_macro = 1;
    flag_conv = 0; flag_ground = false; flag_diis = false;

    %fprintf('\nDIIS Cycle - %d\n', it_macro)
    while it_micro < maxit && flag_conv == 0

        tic
       
        %if mod(it_micro,diis_size) == 0 && it_micro > 1
        if it_micro > diis_size
            flag_diis = true;
            it_macro = it_macro + 1;
            %fprintf('\nDIIS Cycle - %d\n',it_macro)
        end
        
        for j = idx_not_converged
            
            % get current L1 and L2
            l1a = reshape(LAMBDA(sys.posv{1},j),sys.size{1});
            l1b = reshape(LAMBDA(sys.posv{2},j),sys.size{2});
            l2a = reshape(LAMBDA(sys.posv{3},j),sys.size{3});
            l2b = reshape(LAMBDA(sys.posv{4},j),sys.size{4});
            l2c = reshape(LAMBDA(sys.posv{5},j),sys.size{5});

            % build LH diagrammatically with MP denominator separated out
            flag_jacobi = true;
            X1A = build_LH_1A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
            X1B = build_LH_1B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
            X2A = build_LH_2A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
            X2B = build_LH_2B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
            X2C = build_LH_2C(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);

            % update L1 and L2 by Jacobi
            [l1a, l1b, l2a, l2b, l2c] = update_L(l1a,l1b,l2a,l2b,l2c,X1A,X1B,X2A,X2B,X2C,HBar_t,sys,omega(j),shift);
            LAMBDA(:,j) = cat(1,l1a(:),l1b(:),l2a(:),l2b(:),l2c(:));
	
            % buid LH - omega*L residual measure (use full LH)
            [LH] = build_ucc_LHBar(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_ground);
            LAMBDA_resid = LH - omega(j)*LAMBDA(:,j);

            % get lcc energy - returns Ecorr + omega
            E_lcc(j) = sqrt(sum(LH.^2))./sqrt(sum(LAMBDA(:,j).^2)) + Ecorr;

            % check exit condition
            lccsd_resid(j) = norm(LAMBDA_resid);

            % append trial and residual vectors to diis lists
            LAMBDA_list{j}(:,mod(it_micro,diis_size)+1) = LAMBDA(:,j);
            LAMBDA_resid_list{j}(:,mod(it_micro,diis_size)+1) = LAMBDA_resid;

            % diis extrapolate
            if flag_diis
                  LAMBDA(:,j) = diis_xtrap(LAMBDA_list{j},LAMBDA_resid_list{j});
                  flag_diis = false;
            end

        end
        
        % biorthonormalize to right eigenvectors using full matrices
        %LAMBDA = LAMBDA/(LAMBDA'*Rvec(:,1:nroot));
        [LAMBDA, ~] = biorth(LAMBDA,Rvec(:,1:nroot),'left');

        % print status
        fprintf('\nIter-%d:     Elapsed Time: %4.4f\n', it_micro,toc)
        fprintf('-------------------------------------------------------\n')
        for j = 1:nroot
            fprintf('   Root - %d     e = %4.10f     |r| = %4.10f\n',j,E_lcc(j),lccsd_resid(j));
        end


        idx_not_converged = find(lccsd_resid > tol);
        
        if all(lccsd_resid < tol)
            flag_conv = 1;
        end
        
        it_micro = it_micro + 1;
        
    end
    
    % build final residuals
    lccsd_resid = zeros(1,nroot);
    E_lcc = zeros(1,nroot);
    for j = 1:nroot
        % buid LH - omega*L residual measure (use full LH)
        l1a = reshape(LAMBDA(sys.posv{1},j),sys.size{1});
        l1b = reshape(LAMBDA(sys.posv{2},j),sys.size{2});
        l2a = reshape(LAMBDA(sys.posv{3},j),sys.size{3});
        l2b = reshape(LAMBDA(sys.posv{4},j),sys.size{4});
        l2c = reshape(LAMBDA(sys.posv{5},j),sys.size{5});
                
        cc_t.l1a{j+1} = l1a;   cc_t.l1b{j+1} = l1b;
        cc_t.l2a{j+1} = l2a; cc_t.l2b{j+1} = l2b;  cc_t.l2c{j+1} = l2c;
        
        [LH] = build_ucc_LHBar(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_ground);  
        lccsd_resid(j) = norm(LH - omega(j)*LAMBDA(:,j));
        E_lcc(j) = sqrt(sum(LH.^2))./sqrt(sum(LAMBDA(:,j).^2)) + Ecorr;
    end
    
    if flag_conv == 1
        fprintf('\nLeft-EOMUCCSD successfully converged in %d iterations (%4.2f seconds)\n',it_micro,toc(tic_Start));
    else
        fprintf('\nLeft-EOMUCCSD failed to converged in %d iterations\n',maxit)
    end
         
    biorthmat = transpose(conj(LAMBDA))*Rvec(:,1:nroot);
    fprintf('|LR - 1| = %4.12f\n',norm(biorthmat-eye(nroot)))
    
end
