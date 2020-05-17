function [R, omega, r0, res] = eomccsd(HBar,t1,t2,sys,opts)

    fprintf('\n==============++Entering EOM-CCSD Routine++=================\n')
    
    FM = sys.FM; VM = sys.VM; occ = sys.occ; unocc = sys.unocc; 
    nroot = opts.nroot;

    % CCSD dimensions
    [Nocc,Nunocc] = size(HBar{1}{1,2});
    Nov = Nocc*Nunocc; HBar_dim = Nov + Nov^2;
    
    % get Epstein-Nesbet CCSD HBar diagonal
    [D, ~, ~] = HBar_CCSD_diagonal(HBar,t1,t2,sys);
    
    % Matrix-vector product function
    HRmat = @(x) HR_matmat(x,HBar);

    % use CIS guess as initial basis vectors
    if strcmp(opts.init_guess,'cis')
        [omega, B_cis] = get_cis_guess(FM,VM,occ,unocc,nroot,opts.nvec_per_root);
        fprintf('Initial CIS energies:\n')
        for i = 1:length(omega)
            fprintf('      E%d = %4.8f\n',i,omega(i))
        end
        B0 = cat(1,B_cis,zeros(Nov^2,size(B_cis,2)));
        opts.init_guess = 'custom';
        
    % use diagonal basis vector guess (identity)
    else
        B0 = [];
        opts.init_guess = 'diagonal';
    end

    %[ R, L_R, omega, res, it, flag_conv] = davidson_fcn(HRmat, D, HBar_dim, nroot, B0, 'right', opts);
    [R, omega, res, flag_conv] = davidson_fcn_v2(HRmat,D,nroot,B0,opts);

    omega = real(omega);
    
    % calculate r0 for each excited state
    r0 = zeros(nroot,1);
    for i = 1:nroot
        r1 = reshape(R(1:Nov,i),[Nunocc,Nocc]);
        r2 = reshape(R(Nov+1:end,i),[Nunocc,Nunocc,Nocc,Nocc]);
        r0(i) = 1/omega(i)*(einsum(HBar{1}{1,2},r1,'ia,ai->') +...
                            0.25*einsum(HBar{2}{1,1,2,2},r2,'ijab,abij->'));
    end
end

