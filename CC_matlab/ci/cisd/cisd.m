function [Cvec,omega,res] = cisd(nroot,sys)

    Nunocc = sys.Nunocc; Nocc = sys.Nocc;

    opts.maxit = 100;
    opts.tol = 1e-6;
    opts.nvec_per_root = 1;
    opts.max_nvec_per_root = 5;
    opts.flag_verbose = 1;
    opts.thresh_vec = 1e-3;

    Dai = zeros(Nunocc,Nocc);
    Dabij = zeros(Nunocc,Nunocc,Nocc,Nocc);

    for a = 1:Nunocc
        for i = 1:Nocc
            Dai(a,i) = sys.fvv(a,a)-sys.foo(i,i)+sys.Vvoov(a,i,i,a);
            for b = 1:Nunocc
                for j = 1:Nocc
                    Dabij(a,b,i,j) = sys.fvv(a,a) + sys.fvv(b,b) - sys.foo(i,i) - sys.foo(j,j); ...
%                                      +sys.Vvvvv(a,b,a,b) + sys.Voooo(i,j,i,j)...
%                                      +sys.Vvoov(a,i,i,a) + sys.Vvoov(b,j,j,b) ...
%                                      -sys.Vvoov(b,i,i,b) - sys.Vvoov(a,j,j,a);
                end
            end
        end
    end

    D = cat(1,Dai(:),Dabij(:));

    [~,idx] = sort(D,'ascend');

    % Matrix-vector product function
    HCmat = @(x) HC_CISD_matmat(x,sys);

    opts.init_guess = 'diagonal';
    B0 = eye(sys.doubles_dim); B0 = B0(:,idx(1:nroot));

    [ Cvec, eigval, res, flag_conv] = davidson_update_CISD(HCmat, Dai, Dabij, nroot, B0, opts);
    omega = real(eigval);
end

