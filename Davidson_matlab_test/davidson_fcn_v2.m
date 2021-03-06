function [ R, eigval, resid_norm] = davidson_fcn_v2(Ax, D, nroot, opts)

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
    
    e = zeros(nroot,1); e_old = zeros(nroot,1);
    
    % diagonal guess for trial space
    [~,idx_diag] = sort(D,'ascend');
    Id = eye(vec_dim);
    B = Id(:,idx_diag(1:nroot));
        
    curr_size = nroot;
    for it = 1:maxit

        curr_size_old = curr_size;
        e_old = e;

        fprintf('\nIter-%d    Subspace size - %d\n',it,curr_size)
        fprintf('-------------------------------------------------------\n')
        
        B = mgson(B,thresh_vec);
        
        SIGMA = feval(Ax,B);

        G = B'*SIGMA;

        [alpha,GD] = eig(G);
        [e,idx] = sort(diag(GD),'ascend');

        % Right ritz vectors
        V = B*alpha;

        % sorting eigenpairs
        e = e(1:nroot);
        alpha = alpha(:,idx(1:nroot));  % expansion coefficients
        V = V(:,idx(1:nroot));

        % vectorized residual calculation
        RES = SIGMA*alpha - V*diag(e);
        resid_norm = sqrt(sum(RES.^2,1));
        idx_unc = find(resid_norm > tol);
        
        % printing 
        for j = 1:nroot
                fprintf('   Root - %d     e = %4.10f     |r| = %4.10f\n',j,real(e(j)),resid_norm(j));
        end

        Btemp = [];
        for j = 1:length(idx_unc)
            
            J = idx_unc(j);

            RES(:,J) = update_R(RES(:,J),e(J),D);
            RES(:,J) = RES(:,J)/norm(RES(:,J));
            
            r_orth = ortho_root_vec(RES(:,J),[B,Btemp]);
            if norm(r_orth)/norm(RES(:,J)) > thresh_vec
                r_orth = r_orth/norm(r_orth);
                Btemp = cat(2,Btemp,r_orth);
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
            else
                B = cat(2,B,Btemp);
            end
        end
        
        curr_size = size(B,2);
        
        if curr_size == curr_size_old
            fprintf('\nDavidson algorithm stangated without converging... exiting\n')
            R = V; eigval = e;
            return
        end


    end

    fprintf('\nDavidson failed to converge in %d iterations',it)
    R = V; eigval = e; 
        
end
                    
       
function [q] = ortho_root_vec(q,B)

    for i = 1:size(B,2)
        b = B(:,i)/norm(B(:,i));
        q = q - (b'*q)*b;
    end
    
end
