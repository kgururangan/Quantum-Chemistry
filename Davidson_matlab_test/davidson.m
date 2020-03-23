function [ V, W, e, res, flag_conv] = davidson(A, nroot, flag_eig, opts)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

    nvec_per_root = opts.nvec_per_root;
    max_nvec_per_root = opts.max_nvec_per_root;
    maxit = opts.maxit;
    tol = opts.tol;
    flag_verbose = opts.flag_verbose;
    

    fprintf('Beginning Davidson diagonlization algorithm...\n')
    tic

    mat_dim = size(A,1);
    curr_size = nvec_per_root*nroot;
    max_size = max_nvec_per_root*nroot;
    
    e = zeros(nroot,1);
    D = diag(A);
    
    switch opts.init_guess
        
        case 'diagonal'
            
            [~, idx] = sort(D, 'ascend');
            B = eye(mat_dim); B = B(:,idx(1:curr_size));
            
        case 'random'
            
            B = rand(mat_dim,curr_size);
            
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
            
            
        otherwise
            disp('Please enter a valid initial guess type!')
    end
    
    %%%%%%% SOLVE RIGHT EIGENPROBLEM %%%%%%%%%
    if strcmp(flag_eig, 'right')
        
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

            SIGMA = A*B;

            G = B'*SIGMA;

            [alpha,GD] = eig(G);
            [e,idx] = sort(diag(GD),'ascend');

            % Right ritz vectors
            V = B*alpha;
            % Left ritz vectors
            W = B*ctranspose(inv(alpha));
            %W = B/alpha';

            % sorting eigenpairs
            e = e(1:nroot);
            alpha = alpha(:,idx(1:nroot));  % expansion coefficients
            V = V(:,idx(1:nroot));
            W = W(:,idx(1:nroot));

            Q = zeros(mat_dim,nroot);
            resid_norm = zeros(nroot,1);

            for j = 1:nroot

                r = SIGMA*alpha(:,j) - e(j)*V(:,j);
                resid_norm(j) = norm(r);

                Q(:,j) = r/(e(j)-D(j));

                if flag_verbose == 1
                    fprintf('   Root = %d     e = %4.10f     |r| = %4.10f\n',j,real(e(j)),resid_norm(j));
                end
            end

                res = sum(resid_norm);
                eps = sqrt(sum((e(1:nroot)-e_old(1:nroot)).^2));
                if res < tol && eps < tol
                    flag_conv = 1;
                    fprintf('\nDavidson successfully converged in %d iterations (%4.2f seconds)\n',it, toc);
                    %V = B(:,1:nroot);
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

            LAMBDA = AH*B;

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
                    fprintf('\nDavidson successfully converged in %d iterations (%4.2f seconds)\n',it, toc);
                    %V = B(:,1:nroot);
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
                    
       
