function [ Rvec, omega, eom_res, flag_conv] = davidson_solver3(Ax, HBar_t, cc_t, nroot, B0, sys, opts)

% Block-Davidson algorithm for diagonalizing large sparse (non-Hermitian)
% matrices

% References:
% K Hirao and H Nakatusji, J. Comp. Phys. 45, 246-254 (1982)
% R Parrish, E Hohenstein, and T J Martinez, J. Chem. Theory Comput, 12,
% 3003-3007 (2016)

% Input:
% Ax - matrix-vector product function for matrix A to be diagonalized
% D - diagonal of matrix A
% nroot - number of desired eigenpairs
% B0 - initial search space matrix, used in conjunction with
%      opts.init_guess = 'custom'
% opts - struct with fields
%        opts.maxit = maximum # of iterations (e.g. 100)
%        opts.tol = convergence tolerance for each root (e.g. 1e-5)
%        opts.max_nvec_per_root = maximum number of search vectors for each
%        root in subspace (e.g. 10)
%        opts.flag_verbose = print verbosity, 1 or 0

% Output
% R - right eigenvectors
% eigval - eigenvalues
% resid_norm - residuals for each root
% flag_conv - 1 if converged, 0 if not converged, 2 if stagnated


    max_nvec_per_root = opts.max_nvec_per_root;
    maxit = opts.maxit;
    tol = opts.tol;
    flag_verbose = opts.flag_verbose;
    thresh_vec = opts.thresh_vec;
    num_add = opts.nvec_per_root;
    
    if flag_verbose == 1
        fprintf('Entering Davidson diagonlization routine...\n')
        fprintf('Root-by-Root Solver\n')
        fprintf('Settings:\n')
        fprintf('    nroot - %d    maxiter - %d    max subspace dim - %d\n',nroot,maxit,max_nvec_per_root)
        fprintf('    tolerance - %.3fE-06\n', tol*1e6)
        fprintf('    subspace vector threshold - %.3fE-06\n', thresh_vec*1e6)
    end

    tic_Start = tic;

    max_size = max_nvec_per_root;

    % orthonormalize the initial trial space
    B0 = mgson(B0,1e-15);
    
    %%%%%%% SOLVE RIGHT EIGENPROBLEM %%%%%%%%%

    Rvec = [];
    omega = zeros(1,nroot);
    eom_res = zeros(1,nroot);

    for j = 1:nroot

        fprintf('\nEOM-CCSD Root-%d\n',j)
        fprintf('------------------------------------------------------------------------------------------------\n')

        tic_Root = tic;

        B = ortho_root_vec(B0(:,j),Rvec);

        SIGMA = zeros(sys.doubles_dim,max_nvec_per_root);
    
        it = 0; flag_conv = false;
        while it < maxit 

            tic

            if size(B,1) < size(B,2)
                fprintf('WARNING: Number of search vectors greater than dimension of search vectors... MGS will behave erratically!\n')
            end

            curr_size = size(B,2);

            if it > 1
                if curr_size == curr_size_old
                    fprintf('\nDavidson algorithm stangated without converging... exiting\n')
                    break
                end
            end

            sigma = feval(Ax,B(:,curr_size));
            SIGMA(:,curr_size) = sigma;

            G = B'*SIGMA(:,1:curr_size);
            [alpha,GD] = eig(G);
            [e,idx] = sort(diag(GD),'ascend');
           
            % Right ritz vectors
            V = B*alpha;

            % sorting eigenpairs
            e = e(1:num_add);
            alpha = alpha(:,idx(1:num_add));  % expansion coefficients
            V = V(:,idx(1:num_add));

            % vectorized residual calculation
            q = SIGMA(:,1:curr_size)*alpha - V*diag(e);
            resid_norm = sqrt(sum(q.^2));

            % printing for verbose
            fprintf('   Iter - %d     e = %4.10f     |r| = %4.10f     Ndim = %d     Elapsed Time = %4.4f\n',it,real(e),resid_norm,curr_size,toc);

            % check convergence
            if resid_norm < tol
                flag_conv = true;
                break
            end

            q = update_R(q,HBar_t,sys,e,0.0);

            q = q/norm(q);

            % update subspace with preconditioned residual
            if curr_size < max_size
       
                qorth = ortho_root_vec(q,[B,Rvec]);
                if norm(qorth)/norm(q) > thresh_vec
                    qorth = qorth/norm(qorth);
                    B = cat(2,B,qorth);
                end

             else
                fprintf('\nRestarting and collapsing...\n')
                B = B*alpha;
            end

            curr_size_old = curr_size;
            
            it = it + 1;

        end

        Rvec(:,j) = V;
        omega(j) = e;
        eom_res(j) = resid_norm;

        if flag_conv
            fprintf('\nEOM-CCSD successfully converged root %d in %4.4f seconds\n',j,toc(tic_Root))
            
            lccopts.maxit = 50;
            lccopts.diis_size = 5;
            lccopts.tol = 1e-6; 
            lccopts.shift = 0;
            lccopts.nroot = 1;
            [Lvec,~] = leftcc_solver1(omega(j),Rvec(:,j),HBar_t,cc_t,sys,lccopts);
        else
            fprintf('\nEOM-CCSD failed to converge root %d\n',j)
        end

        fprintf('------------------------------------------------------------------------------------------------\n')

    end

    fprintf('\nEOM-CCSD calculation completed in %4.4f seconds\n',toc(tic_Start))

        
    
    
    
end

