function [C, omega, c0] = cis(sys,opts)

    fprintf('\n==================================++Entering CIS Routine++=============================\n')
    
    nroot = opts.nroot;
    
    D1A = zeros(sys.Nvir_alpha,sys.Nocc_alpha);
    D1B = zeros(sys.Nvir_beta,sys.Nocc_beta);
    
    for a = 1:sys.Nvir_alpha
        for i = 1:sys.Nocc_alpha
            D1A(a,i) = sys.fa_vv(a,a) - sys.fa_oo(i,i) + sys.vA_ovvo(i,a,a,i);
        end
    end
    
    for a = 1:sys.Nvir_beta
        for i = 1:sys.Nocc_beta
            D1B(a,i) = sys.fb_vv(a,a) - sys.fb_oo(i,i) + sys.vC_ovvo(i,a,a,i);
        end
    end
    
    D = cat(1,D1A(:),D1B(:));
    [Dsort, idx] = sort(D,'ascend');
    omega_init = Dsort(1:nroot);
    
    fprintf('Initial diagonal energies:\n')
    for i = 1:length(omega_init)
        fprintf('      E%d = %4.8f\n',i,omega_init(i))
    end
    fprintf('\n')
    
    Id = eye(length(D));
    B0 = Id(:,idx(1:nroot));
    
    % Matrix-vector product function
    HRmat = @(x) HC_ucc_matmat(x,sys);
    
    % Davidson eigensolver
    [C, omega, res, ~] = davidson_cis_solver(HRmat, nroot, B0, sys, opts);
    omega = real(omega);
    
    % calculate c0 for each excited state - measures coupling of ground to
    % excited state. Should be 0 by Brillouin's Thm for converged HF
    % orbitals
    c0 = zeros(1,nroot);
    for i = 1:nroot
        c1a = reshape(C(sys.posv{1},i),sys.size{1});
        c1b = reshape(C(sys.posv{2},i),sys.size{2});

        HC0 = +einsum_kg(sys.fa_ov,c1a,'me,em->')...
              +einsum_kg(sys.fb_ov,c1b,'me,em->');

        c0(i) = HC0/omega(i);       
    end


end
