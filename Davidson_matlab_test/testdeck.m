clear all
clc
close all

nroot = 5;
n = 5000;

test_no = 1;

switch test_no
    
    case 1
        % % Test matrix from Hirao, Nakasatuji (non-Hermitian)
        A = zeros(n);
        k = n/2;
        for i = 1:n
            for j = 1:n
                if j <= k
                    A(i,j) = i*(i==j) - (i-j-k^2);
                else
                    A(i,j) = i*(i==j) + (i-j-k^2);
                end
            end
        end

    case 2
        
        sparsity = 1e-5;
        A = diag([1,1,1,2:n-2])+sparsity*rand(n);
        
    case 3
        
        A = full(delsq(numgrid('N',n)));
        
    case 4
        
        % modified Nesbet matrix (Hermitian)
        % from Wood, Zunger J. Phys. A: Math. Gen. 18 (1985) 1343-1359
        A = zeros(n);
        for i = 1:n
            for j = 1:n
                if i ~= j
                    if i < j
                        A(i,j) = 1;
                    else
                        A(i,j) = 1;
                    end
                else
                    if i <= 5
                        A(i,i) = 1 + 0.1*(i-1);
                    else
                        A(i,i) = 2*i-1;
                    end
                end
            end
        end
        
        
    otherwise
        disp('Enter valid test number!')
        
end



%% Exact dense diagonalization

[V_true,D_true, W_true] = eig(A);
[e_true, idx] = sort(diag(D_true),'ascend');
e_true = e_true(1:nroot);
V_true = V_true(:,idx(1:nroot));
W_true = W_true(:,idx(1:nroot));

%% Chebychev davidson

opts=struct( 'blk', 1,  'polym', 20, 'displ', 2, 'tol', 1e-9, 'vimax', 10, 'vomax', 50);  

[e, V]=bchdav(A, nroot, opts); 

%% Vanilla Davidson method

opts.maxit = 1000;
opts.tol = 1e-4;
opts.nvec_per_root = 1;
opts.max_nvec_per_root = 10;
opts.flag_verbose = 1;
opts.init_guess = 'diagonal';
opts.init_space = 5;

[ V, WR, eR, resR, flag_conv] = davidson(A, nroot, 'right', opts);
[ VL, W, eL, resL, flag_conv] = davidson(A, nroot, 'left', opts);

if resR < resL
    e = eR;
else
    e = eL;
end

LR = ctranspose(W)*V;
scale_f = LR\eye(nroot);
V = V*scale_f;

print_errors(V_true, W_true, e_true, V, W, e);

%% Non-Hermitian Davidson method

opts.maxit = 1000;
opts.tol = 1e-4;
opts.nvec_per_root = 1;
opts.max_nvec_per_root = 10;
opts.flag_verbose = 1;
opts.init_guess = 'random';
opts.init_space = 5;

[ V, W, e, flag_conv] = davidson_nonsym(A, nroot, opts);
print_errors(V_true, W_true, e_true, V, W, e);


%%

[y] = Chebyshev_filter(A(:,1), 2, 0, 1, A)

%%

ubberb = upper_bound_lanczos(A,10)

%%
clear all
clc
close all

n = 1000;

A = rand(n);
%A = 0.5*(A+A');
b = rand(n,1);
xt = A\b;

x0 = rand(n,1);
M = eye(n);
restrt = 50;
max_it = 10000;
tol = 1e-6;

[x, error, iter, flag] = gmres_fcn( A, x0, b, M, restrt, max_it, tol );

norm(x-xt)

% C = cell(n,1);
% Ctest = cell(n,1);
% C{1} = chebfun(t,0);
% C{2} = chebfun(t,1);
% 
% for k = 2:n
%     C{k+1} = 2*t.*C{k} - C{k-1};
%     Ctest{k+1} = chebfun(t,k);
%     res = norm(Ctest{k+1}-C{k+1},2);
%     res
% end

%%
clear all
clc
close all

A = rand(8,5) + 1i*rand(8,5);
B = rand(8,5) + 1i*rand(8,5);
%A'*B

[A,B] = biorth(A,B);

ctranspose(A)*B
