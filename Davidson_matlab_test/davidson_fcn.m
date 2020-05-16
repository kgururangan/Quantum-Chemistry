function [ R, L, eigval, res, it, flag_conv] = davidson_fcn(Ax, D, mat_dim, nroot, B0, flag_eig, opts)
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


    nvec_per_root = opts.nvec_per_root;
    max_nvec_per_root = opts.max_nvec_per_root;
    maxit = opts.maxit;
    tol = opts.tol;
    flag_verbose = opts.flag_verbose;
    thresh_vec = opts.thresh_vec;
    
    if flag_verbose == 1
        fprintf('Beginning Davidson diagonlization algorithm...\n')
    end
    
    tic

    curr_size = nvec_per_root*nroot;
    max_size = max_nvec_per_root*nroot;
    
    e = zeros(nroot,1);
    
    switch opts.init_guess
        
        case 'diagonal'
            
            [~, idx] = sort(D, 'ascend');
            B = eye(mat_dim); B = B(:,idx(1:curr_size));
            
        case 'random'
            
            B = rand(mat_dim,curr_size);
            [B,~,~] = mgson(B);
            
        case 'eigen'
            
            n_init = opts.init_space;
            if n_init < nroot
                fprintf('Warning: initial space must be at least as big as number of roots!\n')
                n_init = 2*nroot;
                fprintf('Changing n_init to %d...\n',n_init);
            end
            Asmall = A(1:n_init, 1:n_init);
            [Vs, Ds] = eig(Asmall);
            [ds_sort, idx] = sort(diag(Ds),'ascend');
            B = zeros(mat_dim, n_init);
            for i = 1:n_init
                B(:,i) = cat(1,Vs(:,idx(i)),zeros(mat_dim-n_init,1));
            end
            
        case 'custom'
            
            B = B0;
            if B'*B ~= eye(size(B,2))
                [B,~,~] = mgson(B);
            end
            
            
        otherwise
            disp('Please enter a valid initial guess type!')
    end
    
    %%%%%%% SOLVE RIGHT EIGENPROBLEM %%%%%%%%%
    
    idx_unc = 1:nroot;
    
    if strcmp(flag_eig, 'right')
        
        it = 0; flag_conv = 0; flag_conv_root = zeros(1,nroot);
        while it < maxit && flag_conv == 0

            if size(B,1) < size(B,2)
                fprintf('WARNING: Number of search vectors greater than dimension of search vectors... MGS will behave erratically!\n')
            end

            [B,~,~] = mgson(B,1e-15);
            %[~,B] = qr(B);

            curr_size = size(B,2);
            
            if it > 1
                if curr_size == curr_size_old
                    fprintf('\nDavidson algorithm stangated without converging... exiting\n')
                    flag_conv = 2;
                    break
                end
            end

            e_old = e;

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
            % Left ritz vectors
            W = B*ctranspose(inv(alpha));

            % sorting eigenpairs
            e = e(idx_unc);
            alpha = alpha(:,idx(idx_unc));  % expansion coefficients
            V = V(:,idx(idx_unc));
            W = W(:,idx(idx_unc));

            Q = zeros(mat_dim,length(idx_unc));
            resid_norm = zeros(length(idx_unc),1);

            Btemp = [];
            for j = 1:length(idx_unc)

                r = SIGMA*alpha(:,j) - e(j)*V(:,j);
                resid_norm(j) = norm(r);
                %Q(:,j) = r/(e(j)-D(idx_unc(j)));
                
                q = r/(e(j)-D(j));
                qorth = ortho_root_vec(q,[B,Btemp]);
                if norm(qorth)/norm(q) > thresh_vec
                    Btemp = cat(2,Btemp,qorth/norm(qorth));
                end

                if flag_verbose == 1
                    fprintf('   Root = %d     e = %4.10f     |r| = %4.10f\n',j,real(e(j)),resid_norm(j));
                end
                
                if resid_norm(j) < tol && sqrt((e(j)-e_old(j))^2) < tol && flag_conv_root(j) ~= 1
                    eigval(idx_unc(j)) = e(j);
                    R(:,idx_unc(j)) = V(:,j);
                    L(:,idx_unc(j)) = W(:,j);
                    flag_conv_root(idx_unc(j)) = 1;
                    %idx_unc(j) = [];
                end
                
            end

            res = sum(resid_norm);
            eps = sqrt(sum((e(1:nroot)-e_old(1:nroot)).^2));
            if res < tol && eps < tol
                flag_conv = 1;
                if flag_verbose
                    fprintf('\nDavidson successfully converged in %d iterations (%4.2f seconds)\n',it, toc);
                end
            else
                if curr_size >= max_size
                    if flag_verbose == 1
                        fprintf('\nRestarting and collapsing...\n')
                    end
                    B = B*alpha;
                    e = e_old;
                
                else
                    B = cat(2,B,Btemp);
                    curr_size_old = curr_size;
                end
            end

            it = it + 1;

        end
        
        if flag_conv == 0
            fprintf('\nDavidson failed to converge in %d iterations',it)
            R = V; L = W; eigval = e; 
        end
        
    end
    
    %%%%%%% SOLVE LEFT EIGENPROBLEM %%%%%%%%%
    if strcmp(flag_eig, 'left')
        
        AH = ctranspose(A);
        DH = diag(AH);
    
        it = 0; flag_conv = 0;
        while it < maxit && flag_conv == 0

            if size(B,1) < size(B,2)
                fprintf('WARNING: Number of search vectors greater than dimension of search vectors... MGS will behave erratically!\n')
            end

            [B,~] = mgson(B);

            curr_size = size(B,2);

            e_old = e;

            if flag_verbose == 1
                fprintf('\nIter-%d    Subspace size-%d\n',it,curr_size)
                fprintf('----------------------------------\n')
            end

            LAMBDA = feval(Ax,B);

            G = ctranspose(LAMBDA)*B;

            [~,GD,alpha] = eig(G);
            [e,idx] = sort(diag(GD),'ascend');

            % Left Ritz vector
            W = B*alpha;

            % Right Ritz vector
            V = B*ctranspose(inv(alpha));

            % sorting eigenpairs
            e = e(1:nroot);
            alpha = alpha(:,idx(1:nroot));  % expansion coefficients

            % Left Ritz vector
            W = B*alpha;

            Q = zeros(mat_dim,nroot);
            resid_norm = zeros(nroot,1);

            for j = 1:nroot

                r = LAMBDA*alpha(:,j) - e(j)*W(:,j);
                resid_norm(j) = norm(r);
                

                Q(:,j) = r/(e(j)-DH(j));

                if flag_verbose == 1
                    fprintf('   Root = %d     e = %4.10f     |r| = %4.10f\n',j,real(e(j)),resid_norm(j));
                end
            end

                res = sum(resid_norm);
                eps = sqrt(sum((e(1:nroot)-e_old(1:nroot)).^2));
                if res < tol && eps < tol
                    flag_conv = 1;
                    if flag_verbose
                        fprintf('\nDavidson successfully converged in %d iterations (%4.2f seconds)\n',it, toc);
                    end
                else
                    if curr_size >= max_size
                        if flag_verbose == 1
                            fprintf('\nRestarting and collapsing...\n')
                        end
                        B = B*alpha;
                        e = e_old;
                    else
                        B = cat(2,B,Q);
                    end
                end

                it = it + 1;

        end
    
    end
    
    
end
                    
       
function [q] = ortho_root_vec(q,B)

    for i = 1:size(B,2)
        b = B(:,i)/norm(B(:,i));
        q = q - (b'*q)*b;
    end
    
end
