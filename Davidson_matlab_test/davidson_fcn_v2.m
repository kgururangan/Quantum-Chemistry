function [ R, eigval, resid_norm, flag_conv] = davidson_fcn_v2(Ax, D, nroot, B0, opts)
% Block-Davidson algorithm for diagonalizing large sparse (non-Hermitian)
% matrices
% K Hirao, H Nakatusji. J Comp Phys. 45, 246-254 (1982)

% Input:
% A - matrix to be diagonalized (or function for matrix-vector product)
% nroot - number of desired eigenpairs
% B0 - initial search space matrix, used in conjunction with
%      opts.init_guess = 'custom'
% flag_eig = 'right' or 'left' eigenvectors
% opts - struct with fields
%        opts.maxit = 2000;
%        opts.tol = 1e-4;
%        opts.nvec_per_root = 1;
%        opts.max_nvec_per_root = 10;
%        opts.flag_verbose = 1;
%        opts.init_guess = 'diagonal';
%        opts.init_space = 5;

% Output
% V - right eigenvectors
% W - left eigenvectors
% e - eigenvalues
% res - residual
% it - iterations
% flag_conv - 1 if converged, 0 if not converged


    max_nvec_per_root = opts.max_nvec_per_root;
    maxit = opts.maxit;
    tol = opts.tol;
    flag_verbose = opts.flag_verbose;
    thresh_vec = opts.thresh_vec;
    
    if flag_verbose == 1
        fprintf('Beginning Davidson diagonlization algorithm...\n')
    end
    
    tic

    max_size = max_nvec_per_root*nroot;
    
    e = zeros(nroot,1);
    
    % orthonormalize the initial trial space
    B = mgson(B0,1e-15);
    
    %%%%%%% SOLVE RIGHT EIGENPROBLEM %%%%%%%%%
    
        
    it = 0; flag_conv = 0; 
    while it < maxit && flag_conv == 0

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

        if flag_verbose == 1
            fprintf('\nIter-%d    Subspace size - %d\n',it,curr_size)
            fprintf('-------------------------------------------------------\n')
        end

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
        
        RES = zeros(size(SIGMA,1),nroot);
        resid_norm = zeros(1,nroot);
        idx_unc = []; ct = 1;
        for j = 1:nroot
            RES(:,j) = SIGMA*alpha(:,j) - e(j)*V(:,j);
            resid_norm(j) = norm(RES(:,j));
            if resid_norm(j) > tol
                idx_unc(ct) = j;
                ct = ct + 1;
            end
        end
        
        % printing for verbose 
        if flag_verbose == 1
            for j = 1:nroot
                    fprintf('   Root = %d     e = %4.10f     |r| = %4.10f\n',j,real(e(j)),resid_norm(j));
            end
        end

        Btemp = [];
        for j = 1:length(idx_unc)
            
            J = idx_unc(j);
            %J = j;

            q = RES(:,J)/(e(J)-D(J));
            q = q/norm(q);
            
            qorth = ortho_root_vec(q,[B,Btemp]);
            if norm(qorth)/norm(q) > thresh_vec
                qorth = qorth/norm(qorth);
                Btemp = cat(2,Btemp,qorth);
            end

        end

        if all(resid_norm < tol) 
            flag_conv = 1;
            R = V; eigval = e;
            if flag_verbose
                fprintf('\nDavidson successfully converged in %d iterations (%4.2f seconds)\n',it, toc);
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
                    
       
function [q] = ortho_root_vec(q,B)

    for i = 1:size(B,2)
        b = B(:,i)/norm(B(:,i));
        q = q - (b'*q)*b;
    end
    
end
