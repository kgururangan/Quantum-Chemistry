function [ R, eigval, resid_norm, flag_conv] = davidson_update_R(Ax, Dai, Dabij, nroot, B0, opts)

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


    max_nvec_per_root = opts.max_nvec_per_root;
    maxit = opts.maxit;
    tol = opts.tol;
    flag_verbose = opts.flag_verbose;
    thresh_vec = opts.thresh_vec;
    
    if flag_verbose == 1
        fprintf('Beginning Davidson diagonlization algorithm...\n')
    end

    [Nunocc, Nocc] = size(Dai); 
    Nov = Nunocc*Nocc;
    
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
        idx_not_converged = find(resid_norm > tol);
        
        % printing for verbose 
        if flag_verbose == 1
            for j = 1:nroot
                    fprintf('   Root - %d     e = %4.10f     |r| = %4.10f\n',j,real(e(j)),resid_norm(j));
            end
        end

        Btemp = [];
        for j = 1:length(idx_not_converged)
            
            J = idx_not_converged(j);

            q = RES(:,J);

            q = update_R(reshape(q(1:Nov),Nunocc,Nocc),reshape(q(Nov+1:end),Nunocc,Nunocc,Nocc,Nocc),...
                         e(J),Dai,Dabij);

            q = q/norm(q);
            
            qorth = ortho_root_vec(q,[B,Btemp]);
            if norm(qorth)/norm(q) > thresh_vec
                qorth = qorth/norm(qorth);
                Btemp = cat(2,Btemp,qorth);
            end

        end

        if all(resid_norm <= tol) 
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

function [R] = update_R(r1,r2,e,D1,D2)

    [Nunocc,Nocc] = size(r1);
    
    % don't explicitly antisymmetrize! Seems to mess things up in LCC
    % and even slows down convergence!

    for a = 1:Nunocc
        for i = 1:Nocc
            r1(a,i) = r1(a,i)/(e - D1(a,i));
            for b = 1:Nunocc
                for j = 1:Nocc
                    r2(a,b,i,j) = r2(a,b,i,j)/(e-D2(a,b,i,j));
%                     r2(b,a,i,j) = -r2(a,b,i,j);
%                     r2(a,b,j,i) = -r2(a,b,i,j);
%                     r2(b,a,j,i) = r2(a,b,i,j);
                end
            end
        end
    end



    R = cat(1, r1(:),r2(:));

end
                   
       
function [q] = ortho_root_vec(q,B)

    for i = 1:size(B,2)
        b = B(:,i)/norm(B(:,i));
        q = q - (b'*q)*b;
    end
    
end

function [Q2, Q, R] = mgson(X, varargin)
% Modified Gram-Schmidt orthonormalization (numerical stable version of Gram-Schmidt algorithm) 
% which produces the same result as [Q,R]=qr(X,0)
% Written by Mo Chen (sth4nth@gmail.com).

if isempty(varargin)
    thresh_vec = 1e-15;
else
    thresh_vec = varargin{1};
end

[d,n] = size(X);
m = min(d,n);
R = zeros(m,n);
Q = zeros(d,m);
Q2 = [];
for i = 1:m
    v = X(:,i);
    for j = 1:i-1
        R(j,i) = Q(:,j)'*v;
        v = v-R(j,i)*Q(:,j);
    end
    R(i,i) = norm(v);
    Q(:,i) = v/R(i,i);
    if norm(v) > thresh_vec
        Q2(:,i) = v/R(i,i);
    end
end
R(:,m+1:n) = Q'*X(:,m+1:n);

end
