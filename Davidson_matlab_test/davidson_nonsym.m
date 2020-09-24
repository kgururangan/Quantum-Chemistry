function [ V, W, e, flag_conv] = davidson_nonsym(A, nroot, opts)
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
    AH = ctranspose(A);
    DH = diag(AH);
    
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
    
    BR = B;
    BL = B;
    
    it = 0; flag_conv = 0;
    while it < maxit && flag_conv == 0
        
        if size(BR,1) < size(BR,2)
            fprintf('WARNING: Number of search vectors greater than dimension of search vectors... MGS will behave erratically!\n')
        end
        
        [BR,BL] = biorth(BR,BL);
        
        curr_size = size(BR,2);
        
        e_old = e;
        
        if flag_verbose == 1
            fprintf('\nIter-%d    Subspace size-%d\n',it,curr_size)
            fprintf('----------------------------------\n')
        end
            
        SIGMA = A*BR; 
        LAMBDA = AH*BL;
        
        G = ctranspose(BL)*SIGMA;
        
        [alphaR,GD,alphaL] = eig(G);
        [e,idx] = sort(diag(GD),'ascend');
        
        % sorting eigenpairs
        e = e(1:nroot);
        alphaR = alphaR(:,idx(1:nroot));  % expansion coefficients
        alphaL = alphaL(:,idx(1:nroot));  % expansion coefficients
        
        % Right and left Ritz vectors
        V = BR*alphaR;
        W = BL*alphaL;
        
        QR = zeros(mat_dim,nroot);
        QL = zeros(mat_dim,nroot);
        
        for j = 1:nroot
            
            rR = SIGMA*alphaR(:,j) - e(j)*V(:,j);    normR = norm(rR);
            rL = LAMBDA*alphaL(:,j) - e(j)*W(:,j);   normL = norm(rL);
            
            QR(:,j) = rR/(e(j)-D(j));
            QL(:,j) = rL/(e(j)-DH(j));
            
            if flag_verbose == 1
                fprintf('   Root = %d     e = %4.10f     |R| = %4.10f     |L| = %4.10f\n',j,real(e(j)),normR, normL);
            end
        end
            
            eps = sqrt(sum((e(1:nroot)-e_old(1:nroot)).^2));
            if normR < tol && normL < tol && eps < tol
                flag_conv = 1;
                fprintf('\nDavidson successfully converged in %d iterations (%4.2f seconds)\n',it, toc);
            else
                if curr_size >= max_size
                    if flag_verbose == 1
                        fprintf('\nRestarting and collapsing...\n')
                    end
                    BR = BR*alphaR;
                    BL = BL*alphaL;
                    e = e_old;
                else
                    BR = cat(2,BR,QR);
                    BL = cat(2,BL,QL);
                end
            end
     
            it = it + 1;
            
    end
    
    
end
                    
       
