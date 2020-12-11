function [Lvec,resid,cc_t] = leftcc_solver1(omega,Rvec,HBar_t,cc_t,sys,opts)
% we're using abij ordering for L1 and L2 now!

    fprintf('\n==================================++Entering Left-EOM-UCCSD Routine++===========================\n')
    
    tic_Start = tic;

    tol = opts.tol;
    diis_size = opts.diis_size;
    maxit = opts.maxit;
    shift = opts.excited_shift;
    nroot = opts.nroot;

    %Ecorr = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);

    Lvec = zeros(sys.doubles_dim,nroot);
    resid = zeros(1,nroot); Imat = eye(nroot);

    for j = 1:nroot

        fprintf('\n-------------------- LCCSD - Root %d --------------------\n',j)

        lambda_list = zeros(sys.doubles_dim,diis_size);
        lambda_resid_list = zeros(sys.doubles_dim,diis_size);

        lambda = Rvec(:,j);
        %LAMBDA = rand(sys.doubles_dim,1);

        it_micro = 0; it_macro = 1;
        flag_conv = false; flag_ground = false; 

        tic_Root = tic;

        %fprintf('DIIS Cycle - %d',it_macro)
        while it_micro < maxit
            
            tic

            % get current L1 and L2
            l1a = reshape(lambda(sys.posv{1}),sys.size{1});
            l1b = reshape(lambda(sys.posv{2}),sys.size{2});
            l2a = reshape(lambda(sys.posv{3}),sys.size{3});
            l2b = reshape(lambda(sys.posv{4}),sys.size{4});
            l2c = reshape(lambda(sys.posv{5}),sys.size{5});

            % build LH diagrammatically with MP denominator separated out
            flag_jacobi = true;
            X1A = build_LH_1A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
            X1B = build_LH_1B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
            X2A = build_LH_2A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
            X2B = build_LH_2B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
            X2C = build_LH_2C(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);

            % update L1 and L2 by Jacobi
            [l1a, l1b, l2a, l2b, l2c] = update_L(X1A,X1B,X2A,X2B,X2C,HBar_t,sys,omega(j),shift);
            lambda = cat(1,l1a(:),l1b(:),l2a(:),l2b(:),l2c(:));

            % biorthogonalize to R
            iid = 1:nroot;
            iid(j) = [];    
            
            % REMOVING THIS LINE SEEMED TO HELP RESULTS OF CR-EOMCC(2,3)
            % MATCH UP WITH JUN'S (ON C1 H2O STRETECHED)...
            %lambda = lambda./(lambda'*Rvec(:,j));
            lambda = ortho_root_vec(lambda,Rvec(:,iid));

            
%             LAMBDA = LAMBDA./(LAMBDA'*Rvec(:,j));

	
            % buid LH - omega*L residual measure (use full LH)
            [LH] = build_ucc_LHBar(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_ground);
            LAMBDA_resid = LH - omega(j)*lambda;

            % get lcc energy - returns Ecorr + omega
            E_lcc = sqrt(sum(LH.^2))./sqrt(sum(lambda.^2)); % + Ecorr

            % calculate residual
            lccsd_resid = norm(LAMBDA_resid);

            fprintf('\n    Iter-%d     Residuum = %4.12f      e = %4.12f      Elapsed Time = %4.2f s',it_micro,lccsd_resid,E_lcc,toc);

            % check exit condition
            if lccsd_resid < tol
                flag_conv = true;
                break;
            end
                
            % append trial and residual vectors to diis lists
            lambda_list(:,mod(it_micro,diis_size)+1) = lambda;
            lambda_resid_list(:,mod(it_micro,diis_size)+1) = LAMBDA_resid;

            % diis extrapolate
            %if mod(it_micro,diis_size) == 0 && it_micro > 1
            if it_micro > diis_size
               it_macro = it_macro + 1;
               %fprintf('\nDIIS Cycle - %d',it_macro)
               lambda = diis_xtrap(lambda_list,lambda_resid_list);
            end        

            it_micro = it_micro + 1;

        end

        Lvec(:,j) = lambda;
        resid(j) = lccsd_resid;

        biorthvec = Lvec(:,j)'*Rvec(:,1:nroot);
        fprintf('\nBiorthogonality measure: |LR - 1| = %4.12f\n',norm(biorthvec-Imat(j,:)))

        if flag_conv
            fprintf('\nLeft-CC for root %d successfully converged in %d iterations (%4.4f seconds)\n',j,it_micro,toc(tic_Root))
        else
            fprintf('\nLeft-CC failed to converge for root %d\n',j)
        end
        
        cc_t.l1a{j+1} = reshape(Lvec(sys.posv{1},j),sys.size{1});
        cc_t.l1b{j+1} = reshape(Lvec(sys.posv{2},j),sys.size{2});
        cc_t.l2a{j+1} = reshape(Lvec(sys.posv{3},j),sys.size{3});
        cc_t.l2b{j+1} = reshape(Lvec(sys.posv{4},j),sys.size{4});
        cc_t.l2c{j+1} = reshape(Lvec(sys.posv{5},j),sys.size{5});
        
        fprintf('---------------------------------------------------------\n')
    end

    LR = einsum_kg(Lvec,Rvec,'pr,ps->rs');
    
    fprintf('\nLeft-CC completed in %4.4fs\n|LR| = %4.12f\n',toc(tic_Start),norm(LR))
    
         
    
end
