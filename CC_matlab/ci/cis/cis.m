function [Cvec,omega,res] = cis(nroot,sys)

    Nunocc = sys.Nunocc/2; Nocc = sys.Nocc/2; Norb = sys_cc.Norb/2;

    opts.maxit = 100;
    opts.tol = 1e-6;
    opts.nvec_per_root = 1;
    opts.max_nvec_per_root = 5;
    opts.flag_verbose = 1;
    opts.thresh_vec = 1e-5;

    Dai = zeros(Nunocc,Nocc);

    for a = Nocc+1:Norb
        for i = 1:Nocc
            Dai(a,i) = sys.e1int(a,a)-sys.e1int(i,i)+sys.e2int(a,i,i,a);
        end
    end

    D = Dai(:);

    [~,idx] = sort(D,'ascend');

    % Matrix-vector product function
    HCmat = @(x) HC_matmat(x,sys);

    opts.init_guess = 'diagonal';
    B0 = eye(sys.singles_dim); B0 = B0(:,idx(1:nroot));

    [ Cvec, eigval, res, flag_conv] = davidson_update_C(HCmat, Dai, nroot, B0, opts);
    omega = real(eigval);
end

