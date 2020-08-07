function [R, omega, r0, res] = eomccsd(HBar,t1,t2,sys,sys_ucc,opts)

    fprintf('\n==================================++Entering EOM-CCSD Routine++=============================\n')
    
    nroot = opts.nroot;

    Nocc = sys.Nocc; Nunocc = sys.Nunocc;
    
    % get Epstein-Nesbet CCSD HBar diagonal
    %[D, ~, ~] = HBar_CCSD_SD_diagonal(HBar,t1,t2,sys);

    Dai = zeros(Nunocc,Nocc);
    Dabij = zeros(Nunocc,Nunocc,Nocc,Nocc);

    % approximate Epstein-Nesbet denominator by one-body HBar
    % <ia|H1|ia> and <ijab|H1|ijab>
    for a = 1:Nunocc
        for i = 1:Nocc
            [D1,D2,D3] = HBar_CCSD_S_diagonal(a,i,HBar);
            Dai(a,i) = D1;
            for b = 1:Nunocc
                for j = 1:Nocc
                    [D1,D2,D3] = HBar_CCSD_D_diagonal(a,b,i,j,HBar);
                    Dabij(a,b,i,j) = D1;
                end
            end
        end
    end

    D = cat(1, Dai(:), Dabij(:));

%     
    % Matrix-vector product function
    HRmat = @(x) HR_matmat(x,HBar);

    % use CIS guess as initial basis vectors
    if strcmp(opts.init_guess,'cis')
% 
%         oa = sys.occ_list(1:2:end); ob = sys.occ_list(2:2:end);
%         ua = sys.unocc_list(1:2:end); ub = sys.unocc_list(2:2:end);
%         ua = ua - sys.occ_list(end); ub = ub - sys.occ_list(end);
% 
%         % use CIS guess as initial basis vectors
%         [omega, C1] = cis_spinadapt(nroot,sys_ucc,1);
% 
%         % closed-shell singlet spin-adaption sets r1a = r1b
%         B1A = C1;
%         B1B = C1;
% 
%         B2A = zeros(length(sys_ucc.posv{3}),size(B1A,2));
%         B2B = zeros(length(sys_ucc.posv{4}),size(B1A,2));
%         B2C = zeros(length(sys_ucc.posv{5}),size(B1A,2));
% 
%         B0 = zeros(sys.doubles_dim,size(B1A,2));
%         for i = 1:size(B1A,2)
%             b1a = reshape(B1A(:,i),sys_ucc.size{1});
%             b1b = reshape(B1B(:,i),sys_ucc.size{2});
%             b2a = reshape(B2A(:,i),sys_ucc.size{3});
%             b2b = reshape(B2B(:,i),sys_ucc.size{4});
%             b2c = reshape(B2C(:,i),sys_ucc.size{5});
%             b1 = convert_spinint_to_spinorb({b1a,b1b},sys_ucc);
%             b2 = convert_spinint_to_spinorb({b2a,b2b,b2c},sys_ucc);
%             B0(:,i) = cat(1,b1(:),b2(:));
%         end

        [omega, B_cis] = get_cis_guess(sys,nroot,opts.nvec_per_root);
        B0 = cat(1,B_cis,zeros(sys.Nocc^2*sys.Nunocc^2,size(B_cis,2)));
        opts.init_guess = 'custom';

        fprintf('Initial CIS energies:\n')
        for i = 1:length(omega)
            fprintf('      E%d = %4.8f\n',i,omega(i))
        end
        fprintf('\n')

    elseif strcmp(opts.init_guess,'cisd')

        fprintf('\n\n==================++Entering CISD Subroutine++====================\n')

        [B0,omega,~] = cisd(nroot,sys);

        fprintf('Initial CISD energies:\n')
        for i = 1:length(omega)
            fprintf('      E%d = %4.8f\n',i,omega(i))
        end
        fprintf('\n')

        opts.init_guess = 'custom';
        
    % use diagonal basis vector guess (identity)
    else
        B0 = [];
        opts.init_guess = 'diagonal';
    end

    %[ R, L_R, omega, res, it, flag_conv] = davidson_fcn(HRmat, D, HBar_dim, nroot, B0, 'right', opts);
    %[R, omega, res, flag_conv] = davidson(HRmat,D,nroot,B0,opts);
    [R, omega, res, flag_conv] = davidson_update_R(HRmat,Dai,Dabij,nroot,B0,opts);

    omega = real(omega);
    
    % calculate r0 for each excited state
    r0 = zeros(nroot,1);
    for i = 1:nroot
        r1 = reshape(R(sys.posv1,i),[Nunocc,Nocc]);
        r2 = reshape(R(sys.posv2,i),[Nunocc,Nunocc,Nocc,Nocc]);
        r0(i) = 1/omega(i)*(einsum_kg(HBar{1}{1,2},r1,'ia,ai->') +...
                            0.25*einsum_kg(HBar{2}{1,1,2,2},r2,'ijab,abij->'));
    end
    r0 = real(r0);
end

