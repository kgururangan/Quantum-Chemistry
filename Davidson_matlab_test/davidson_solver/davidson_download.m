function [V, E, FLAG, RES] = davidson_download(A, nroot, iter_count_max, tol, varargin)
%DAVIDSON Davidson iterative method for a few of the eigenvalues and 
%eigenvectors of real symmetric and diagonally dominant matrices
% E = DAVIDSON(A) returns the lowest eigenvalue of a real symmetric and
% diagonally dominant matrix A.
% E = DAVIDSON(A, I) returns the lowest few eigenvalues specified by index
% array I.
% E = DAVIDSON(..., 'largest') returns the largest (few) eigenvalue(s).
% E = DAVIDSON(..., 'rough') uses a less strict convergence criterion.
% [V, E] = DAVIDSON(...) returns the eigenpairs.
% [V, E, FLAG] = DAVIDSON(...) also returns a convergence flag. FLAG is the
% number of eigenvalues that are not converged.

% Reference:
% [1] Davidson, E.R., "The iterative calculation of a few of the lowest 
% eigenvalues and corresponding eigenvectors of large real-symmetric 
% matrices", J. Comput. Phys. 17, 87-94 (1975)

debug = 0;

if isempty(varargin)
    opt_eval = 'smallest';
else
    opt_eval = varagin{1};
end

if strcmp(opt_eval,'largest')
    A = -A;
end


I = 1:nroot;

sz = size(A,1);
if max(I) > sz
    error('Index exceeds matrix size.')
end

[I_sort, idx_sort_I] = sort(I);
if I_sort(end) > nthroot(sz, 3)
    warning('Computing general eigenvalues is inefficient.')
end

nvals = length(I); % number of eigenvalues to compute
E = zeros(nvals, 1); % eigenvalues
V = zeros(sz, nvals); % eigenvectors

% convergence criteria
crit1 = @(q) norm(q) < tol;
crit2 = @(q, alphaM_end) norm(q) < 1e-4 && abs(alphaM_end) < 1e-12;

% max number of iterations for each eigenvalue

if debug
    hist_conv_crit = zeros(nvals, 1);
    hist_iter_count = zeros(nvals, 1);
    hist_qnorm = zeros(iter_count_max, nvals);
    hist_alpha = zeros(iter_count_max, nvals);
end

% subspace size limits for the k-th eigenvalue
sz_subspace_min = @(k) k + 5;
sz_subspace_max = @(k) k + 20;

if sz_subspace_max(I_sort(end)) > sz
    warning('Matrix size is too small. Use built-in dense algorithm instead.')
    [V_all, E_all] = eig(A, 'vector');
    [E_all, idx_sort_E_all] = sort(E_all);
    V_all = V_all(:, idx_sort_E_all);
    
    if strcmp(opt_eval, 'largest')
        E_all = -E_all;
    end
    
    E = E_all(I);
    V = V_all(I);
    
    if nargout == 1
        V = E;
    end
    
    FLAG = sum(isnan(E));
    return
end

% preallocate space for b and Ab
b = zeros(sz, sz_subspace_min(I_sort(end)));
Ab = zeros(sz, sz_subspace_max(I_sort(end)));

RES = zeros(1,length(I));

for ival = 1 : nvals
    
    fprintf('\nSOLVING FOR ROOT %d...\n',ival)
    fprintf('----------------------------------------------------------------\n')
    
    k = I_sort(ival); % seek the k-th eigenvalue
    
    % subspace size limits for the k-th eigenvalue
    sz_min = sz_subspace_min(k);
    sz_max = sz_subspace_max(k);
    
    % guess space
    if ival == 1
        % generate from scratch
        b(:, 1:sz_min) = guess_init(diag(A), sz_min);
        sz_b = sz_min;
    else
        sz_alpha = size(alphaM, 1);
        b(:, 1:sz_alpha) = b(:, 1:sz_alpha) * alphaM;
        sz_b = sz_alpha;
        
        if sz_b < sz_min
            b(:, 1:sz_min) = [b(:, 1:sz_b), guess_init(diag(A), sz_min-sz_b)];
            sz_b = sz_min;
        end
        [b(:, 1:sz_b), ~] = qr(b(:, 1:sz_b), 0);
    end
    
    Ab(:, 1:sz_b) = A * b(:, 1:sz_b);
    
    iter_count = 0;
    while true
                
        if iter_count > iter_count_max
            E(ival) = nan;
            V(:, ival) = nan;
            if debug
                hist_iter_count(ival) = iter_count;
            end
            fprintf('FAILED TO CONVERGE ROOT %d\n',k)
            break
        end
        Asub = b(:, 1:sz_b)' * Ab(:, 1:sz_b);
        %Asub = (Asub + conj(Asub')) / 2; % enforces hermiticity of matrix
        
        [alphaM, lambda] = eig(Asub, 'vector');
        [lambda, idx_sort_lambda] = sort(lambda);
        alphaM = alphaM(:, idx_sort_lambda);

        if sz_b == sz_max
            b(:, 1:sz_min) = b(:, 1:sz_max) * alphaM(:, 1:sz_min);
            sz_b = sz_min;
            [b(:, 1:sz_b), ~] = qr(b(:, 1:sz_b), 0);
            Ab(:, 1:sz_b) = A * b(:, 1:sz_b);
            continue
        end
        
        % residual vector
        q = ( Ab(:, 1:sz_b) - lambda(k)*b(:, 1:sz_b) ) * alphaM(:, k);
        
        fprintf('Iter-%d    |r| = %4.8f    e = %4.8f    Subspace size - %d\n',iter_count,norm(q),lambda(k),sz_b)
        
        %is_conv = false;
        %if strcmp(opt_conv, 'strict')
        is_conv = crit1(q);
        %elseif strcmp(opt_conv, 'rough')
        %    is_conv = crit1(q) || crit2(q, alphaM(end, k));
        %end
        
        if is_conv
            E(ival) = lambda(k);
            RES(ival) = norm(q);
            V(:, ival) = b(:, 1:sz_b) * alphaM(:, k);
            if debug
                hist_conv_crit(ival) = 1*crit1(q) + 2*crit2(q, alphaM(end, k));
                hist_iter_count(ival) = iter_count;
            end
            fprintf('SUCCESSFULLY CONVERGED ROOT %d\n',ival)
            break
        end
        
        d = q ./ (lambda(k) - diag(A));
        if any(~isfinite(d))
            d(isfinite(d)) = 0;
            d(~isfinite(d)) = 1;
        end
        [b_new, ~] = qr([b(:, 1:sz_b), d], 0);
        b(:, sz_b+1) = b_new(:, end);
        Ab(:, sz_b+1) = A*b_new(:, end);
        sz_b = sz_b + 1;
        
        iter_count = iter_count + 1;
        
        if debug
            hist_qnorm(iter_count, ival) = norm(q);
            hist_alpha(iter_count, ival) = abs(alphaM(end, k));
        end
    end
end

if strcmp(opt_eval, 'largest')
    E = -E;
end

[~, idx_revsort_I] = sort(idx_sort_I);
E = E(idx_revsort_I);
V = V(:, idx_revsort_I);

if nargout == 1
    V = E;
end

FLAG = sum(isnan(E));

end

function [ b ] = guess_init(A_diag, sz_sub)
% generate initial guess space based on the the diagonal elements
sz = length(A_diag);
[~, idx_small] = mink(A_diag, sz_sub);
b = 0.001 * (rand(sz, sz_sub) - 0.5);
%b = 0.001 * ones(sz, sz_sub);
b(idx_small + (0:sz_sub-1)'*sz) = 1;
[b, ~] = qr(b, 0);
end

