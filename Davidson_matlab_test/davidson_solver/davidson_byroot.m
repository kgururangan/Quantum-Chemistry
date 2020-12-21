function [ R, eigval, resid_norm] = davidson_byroot(Ax, D, nroot, opts)

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
    
    R = zeros(vec_dim,nroot);
    eigval = zeros(1,nroot);
    
    % diagonal guess for trial space
    [~,idx_diag] = sort(D,'ascend');
    Id = eye(vec_dim); Id(:,idx_diag(1:nroot));
    
    for ival = 1:nroot
        
            fprintf('SOLVING FOR ROOT %d\n',ival)
            fprintf('--------------------------------------------------------------\n')
    
            B = zeros(vec_dim,nvec);
            SIGMA = zeros(vec_dim,nvec);
            B(:,1) = Id(:,ival);

            curr_size = 1;
            for it = 1:maxit

                B0 = mgson(B(:,1:curr_size));
                B(:,1:curr_size) = B0;

                sigma = feval(Ax,B(:,curr_size));
                SIGMA(:,curr_size) = sigma;

                G = B(:,1:curr_size)'*SIGMA(:,1:curr_size);
                [alpha,GD] = eig(G);
                [e,idx] = sort(diag(GD),'ascend');

                % sorting eigenpairs
                e = e(1);
                alpha = alpha(:,idx(1));  % expansion coefficients
                V = B(:,1:curr_size)*alpha;

                % vectorized residual calculation
                RES = SIGMA(:,1:curr_size)*alpha - V*e;
                resid_norm = norm(RES);

                % printing 
                fprintf('   Iter - %d     e = %4.10f     |r| = %4.10f\n',it,real(e),resid_norm);

                RES = update_R(RES,e,D);
                RES = RES/norm(RES);
                r_orth = orthogonalize_root(RES,B(:,1:curr_size));
                %norm(r_orth)
                if norm(r_orth)/norm(RES) > thresh_vec
                    add_B = r_orth/norm(r_orth);
                end

                if resid_norm < tol
                    R(:,ival) = V;  eigval(ival) = e;
                    fprintf('\nSUCCESSFULLY CONVERGED ROOT %d (%4.2f seconds)\n',ival, toc);
                    break
                else
                    if curr_size >= nvec
                        fprintf('\nRestarting and collapsing...\n')
                        B = B*alpha;
                        curr_size = 1;
                    else
                        B(:,curr_size+1) = add_B;
                        curr_size = curr_size +  1;
                    end
                end


            end
            
            fprintf('\nFAILED TO CONVERGE ROOT %d',ival)
    
    end
        
end
