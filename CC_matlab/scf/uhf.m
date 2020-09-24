function [Escf, C, P, eps_mo, Fock_A] = uhf(sys,opts,P_guess)

    Smat = sys.overlap;
    Zmat = sys.e1int;
    VVmat = sys.e2int;
    Vnuc = sys.Vnuc;
    
    diis_size = opts.diis_size;
    level_shift = opts.level_shift;
    maxit = opts.maxit;
    tol = opts.tol;
    
    fprintf('\n==================================++Entering UHF Routine++=================================\n')
    
    Norb = size(Smat,1);    
    Nocc_a = sys.Nocc_a;
    Nocc_b = sys.Nocc_b;

    tic_Start = tic;
    
    % Orthogonalization of basis
    X = orthomat(Smat);
    
    % Initial guess at P 
    if nargin < 3
        Pa = zeros(Norb);
        Pb = zeros(Norb) + 1e-3*rand(Norb);
    else
        Pa = P_guess{1};
        Pb = P_guess{2};
    end
    
%     % Calculate Fock matrix with guess
    Ga = einsum_kg(VVmat,Pa,'prqs,rs->pq') - ...
        0.5*einsum_kg(VVmat,Pa,'prsq,rs->pq');    
    Fock_A = Zmat + Ga;% + level_shift*(Smat-Smat.*P.*Smat);
    
    % Diagonalize Fock and get MO vectors and energies
    [Cp,eps_mo] = eig(X'*Fock_A*X); eps_mo = diag(eps_mo);
    [~,idx] = sort(eps_mo,'ascend'); eps_mo = eps_mo(idx);   
    C = X*Cp(:,idx);
    
    % Build density matrix
    P = 2*einsum_kg(C(:,1:Nocc),C(:,1:Nocc),'pi,qi->pq');

    Escf = 0.5*einsum_kg(P,Zmat+Fock_A,'pq,qp->') + Vnuc;

    % DIIS containers
    F_list = zeros(Norb^2,diis_size); F_resid_list = zeros(Norb^2,diis_size);
    
    % SCF Iterations
    it_micro = 0; it_macro = 0; flag_conv = false;
    while it_micro < maxit
        
        tic
        
        % Calculate Fock matrix 
        Ga = einsum_kg(VVmat,P,'prqs,rs->pq') - ...
            0.5*einsum_kg(VVmat,P,'prsq,rs->pq');

        Fock_A = Zmat + Ga + level_shift*(Smat-Smat*P*Smat);

        % Build DIIS Residual
        diis_r = X*(Fock_A*P*Smat - Smat*P*Fock_A)*X;
        scf_resid = sqrt(mean(diis_r(:).^2));
      
        % append trial and residual vectors to lists
        F_list(:,mod(it_micro,diis_size)+1) = Fock_A(:);
        F_resid_list(:,mod(it_micro,diis_size)+1) = diis_r(:);
         
        % diis extrapolate
        if mod(it_micro,diis_size) == 0 && it_micro > 0
           it_macro = it_macro + 1;
           fprintf('\nDIIS Cycle - %d',it_macro)
           Fock_A = diis_xtrap(F_list,F_resid_list);
           Fock_A = reshape(Fock_A,[Norb,Norb]);
        end    
        
        % Diagonalize Fock and get MO vectors and energies
        [Cp,eps_mo] = eig(X'*Fock_A*X); eps_mo = diag(eps_mo);
        [~,idx] = sort(eps_mo,'ascend'); eps_mo = eps_mo(idx);
        C = X*Cp(:,idx);

        % Compute HF energy and density matrix 
        P = 2*einsum_kg(C(:,1:Nocc),C(:,1:Nocc),'pi,qi->pq');
        Escf = 0.5*einsum_kg(P,Zmat+Fock_A,'qp,pq->') + Vnuc;

        fprintf('\n    Iter-%d     Residuum = %4.12f      Escf = %4.12f      Elapsed Time = %4.2f s',it_micro,scf_resid,Escf,toc);
        
        % Check convergence condition
        if scf_resid <= tol
            flag_conv = true;
            break
        end
        
        it_micro = it_micro + 1;
                

    end

    
    if flag_conv
        fprintf('\nRHF successfully converged in %d iterations (%4.2f seconds)\n',it_micro,toc(tic_Start));
        fprintf('Total Energy = %4.12f Ha\n',Escf);
    else
        fprintf('\nRHF failed to converged in %d iterations\n',maxit)
    end
    
end

