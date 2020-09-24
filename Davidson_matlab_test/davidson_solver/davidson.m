function [ R, eigval, resid_norm] = davidson(Ax, D, nroot, opts)

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
% V - right eigenvectors
% e - eigenvalues
% res - residuals for each root
% flag_conv - 1 if converged, 0 if not converged, 2 if stagnated


    nvec = opts.nvec;
    maxit = opts.maxit;
    tol = opts.tol;
    thresh_vec = opts.thresh_vec;

    fprintf('Beginning Davidson diagonlization algorithm...\n')

    
    tic
    
    
    %%%%%%% SOLVE RIGHT EIGENPROBLEM %%%%%%%%%
       

    max_size = nvec*nroot;
    vec_dim = length(D);
    
    e = zeros(nroot,1); 
    
    % diagonal guess for trial space
    [~,idx_diag] = sort(D,'ascend');
    Id = eye(vec_dim);
    B0 = Id(:,idx_diag(1:nroot));
    
    B = zeros(vec_dim,max_size);
    SIGMA = zeros(vec_dim,max_size);
    B(:,1:nroot) = B0;
        
    curr_size = nroot;
    for it = 1:maxit

        if it == 1
        	n_prev = 1;
        else
            n_prev = curr_size - ct_add + 1;
        end

        fprintf('\nIter-%d    Subspace size - %d\n',it,curr_size)
        fprintf('-------------------------------------------------------\n')
        
        B0 = mgson(B(:,1:curr_size));
        %B0 = gramschmidt(B(:,1:curr_size));
        B(:,1:curr_size) = B0;
        
        sigma = feval(Ax,B(:,n_prev:curr_size));
        SIGMA(:,n_prev:curr_size) = sigma;
        
        G = B(:,1:curr_size)'*SIGMA(:,1:curr_size);
        [alpha,GD] = eig(G);
        [e,idx] = sort(diag(GD),'ascend');

        % Right ritz vectors
        V = B(:,1:curr_size)*alpha;

        % sorting eigenpairs
        e = e(1:nroot);
        alpha = alpha(:,idx(1:nroot));  % expansion coefficients
        V = V(:,idx(1:nroot));

        % vectorized residual calculation
        RES = SIGMA(:,1:curr_size)*alpha - V*diag(e);
        resid_norm = sqrt(sum(RES.^2,1));
        idx_unc = find(resid_norm > tol);
        
        % printing 
        for j = 1:nroot
                fprintf('   Root - %d     e = %4.10f     |r| = %4.10f\n',j,real(e(j)),resid_norm(j));
        end

        add_B = zeros(vec_dim,nroot); ct_add = 0;
        for j = 1:length(idx_unc)
            
            J = idx_unc(j);

            RES(:,J) = update_R(RES(:,J),e(J),D);
            RES(:,J) = RES(:,J)/norm(RES(:,J));
            
            if ct_add > 0
                r_orth = orthogonalize_root(RES(:,J),[B(:,1:curr_size),add_B(:,1:ct_add)]);
            else
                r_orth = orthogonalize_root(RES(:,J),B(:,1:curr_size));
            end
            
            if norm(r_orth)/norm(RES(:,J)) > thresh_vec
                ct_add = ct_add + 1;
                add_B(:,ct_add) = r_orth/norm(r_orth);
            end

        end

        if all(resid_norm <= tol)
            R = V; eigval = e;
            fprintf('\nDavidson successfully converged in %d iterations (%4.2f seconds)\n',it, toc);
            return
        else
            if curr_size >= max_size
                fprintf('\nRestarting and collapsing...\n')
                B = B*alpha;
                curr_size = nroot;
            else
                B(:,curr_size+1:curr_size+ct_add) = add_B(:,1:ct_add);
                curr_size = curr_size + ct_add;
            end
        end


    end

    fprintf('\nDavidson failed to converge in %d iterations',it)
    R = V; eigval = e; 
        
end

