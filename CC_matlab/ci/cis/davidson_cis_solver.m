function [ R, eigval, resid_norm, flag_conv] = davidson_cis_solver(Ax, nroot, B0, sys, opts)

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
    
    if flag_verbose == 1
        fprintf('Entering Davidson diagonlization routine...\n')
        fprintf('Block Davidson Solver\n')
        fprintf('Settings:\n')
        fprintf('    nroot - %d    maxiter - %d    max subspace dim - %d\n',nroot,maxit,max_nvec_per_root*nroot)
        fprintf('    tolerance - %.3fE-06\n', tol*1e6)
        fprintf('    subspace vector threshold - %.3fE-06\n', thresh_vec*1e6)
    end

    tic_Start = tic;

    max_size = max_nvec_per_root*nroot;

    % orthonormalize the initial trial space
    B = mgson(B0,1e-15);
    
    %%%%%%% SOLVE RIGHT EIGENPROBLEM %%%%%%%%%
 
    SIGMA = zeros(sys.singles_dim,max_size);
        
    it = 0; flag_conv = 0; 
    while it < maxit && flag_conv == 0

    tic

        if size(B,1) < size(B,2)
            fprintf('WARNING: Number of search vectors greater than dimension of search vectors... MGS will behave erratically!\n')
        end

        curr_size = size(B,2);

        if it > 1
            if curr_size == curr_size_old
                fprintf('\nDavidson algorithm stangated without converging... exiting\n')
                flag_conv = 2;
                break
            end
        end
        
        % orthonormalize search space using modified gram-schmidt
        B = mgson(B,thresh_vec);

        % evaluate HC sigma product for new search vectors (the last nroot of them)
        sigma = feval(Ax,B(:,max(curr_size-nroot+1,1):curr_size));
        SIGMA(:,max(curr_size-nroot+1,1):curr_size) = sigma;
        G = B'*SIGMA(:,1:curr_size);
        
        % solve small projection subspace eigenvalue problem
        [alpha,GD] = eig(G);
        [e,idx] = sort(diag(GD),'ascend');

        % Right ritz vectors
        V = B*alpha;

        % sorting eigenpairs
        e = e(1:nroot);
        alpha = alpha(:,idx(1:nroot));  % expansion coefficients
        V = V(:,idx(1:nroot));
        
        % vectorized residual calculation
        RES = SIGMA(:,1:curr_size)*alpha - V*diag(e);
        resid_norm = sqrt(sum(RES.^2,1));
        idx_not_converged = find(resid_norm > tol);

        Btemp = [];
        for j = 1:length(idx_not_converged)
        %for J = 1:nroot
            
            J = idx_not_converged(j);

            q = RES(:,J);

            q = update_C(q,sys,e(J),0.0);

            q = q/norm(q);
            
            qorth = ortho_root_vec(q,[B,Btemp]);
            if norm(qorth)/norm(q) > thresh_vec
                qorth = qorth/norm(qorth);
                Btemp = cat(2,Btemp,qorth);
            end

        end


        % printing for verbose
        if flag_verbose == 1
            fprintf('\nIter-%d    Subspace size - %d    Elapsed Time - %4.4fs\n',it,curr_size,toc)
            fprintf('-------------------------------------------------------\n')
            for j = 1:nroot
                    fprintf('   Root - %d     e = %4.10f     |r| = %4.10f\n',j,real(e(j)),resid_norm(j));
            end
        end

        % check convergence
        if all(resid_norm <= tol) 
            flag_conv = 1;
            R = V; eigval = e;
            if flag_verbose
                fprintf('\nDavidson successfully converged in %d iterations (%4.2f seconds)\n',it, toc(tic_Start));
            end
        else
            if curr_size >= max_size
                if flag_verbose == 1
                    fprintf('\nRestarting and collapsing...\n')
                end
                B = B*alpha;

            else
                B = cat(2,B,Btemp);
                curr_size_old = curr_size;
            end
        end

        it = it + 1;

    end

    if flag_conv == 0
        fprintf('\nDavidson failed to converge in %d iterations',it)
        R = V; eigval = e; 
    end
        
    
    
    
end
