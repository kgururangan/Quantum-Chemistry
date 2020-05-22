function [Fock, C, P, EHF] = scf(Zmat, VVmat, Smat)

    addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab/utils'));


    Norb = size(Zmat,1);
    Nelec = Norb;
    Nocc0 = Nelec/2;

    % orthonormalization matrix wrt to overlap 
    X = orthomat(Smat,'canonical');

    % SCF procedure  
    P = zeros(Norb); 
    G = einsum(VVmat,P,'pqrs,rs->pq') - 0.5*einsum(VVmat,P,'prsq,rs->pq');
    Fock = Zmat + G;
    
    [C, D_mo] = eig(X'*Fock*X);
    [eps_mo,idx] = sort(diag(D_mo),'ascend');
    C = X*C(:,idx);
    Fock = C'*Fock*C;
    P0 = 2*einsum(C(:,1:Nocc0),C(:,1:Nocc0),'pi,qi->pq');
    P = P0;
    
    tol = 1e-9; it = 1; maxit = 1e3; flag_exit = 0;
    while flag_exit == 0 
        
        % calculate G
        G = einsum(VVmat,P,'pqrs,rs->pq') - 0.5*einsum(VVmat,P,'prsq,rs->pq');
        Fock = Zmat + G;
    
        % diagonalize Fock matrix
        [C, D_mo] = eig(X'*Fock*X);
        [eps_mo, idx] = sort(diag(D_mo),'ascend');
        
        % obtain molecular orbital transformation coeffs, order by energy
        C = X*C(:,idx);
        Fock = C'*Fock*C;
        
        % calculate density matrix
        P = 2*einsum(C(:,1:Nocc0),C(:,1:Nocc0),'pi,qi->pq');
        
        % check exit condition
        dp_rms = sqrt(mean((P(:)-P0(:)).^2));
        if dp_rms < tol
            flag_exit = 1;
        elseif it > maxit
            flag_exit = 2;
        else
            P0 = P;
            it = it + 1;
        end

    end

    % Matrix elements in molecular orbital (Fock) basis
    EHF = 0.5*einsum(P,Zmat+Fock,'qp,pq->');
    
    if flag_exit == 1
        fprintf('SCF successfully converged in %d interations\n',it);
    else
        fprintf('SCF failed to converge in %d iterations\n',it);
    end
end

