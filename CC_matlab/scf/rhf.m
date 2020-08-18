function [Escf, C, P, eps_mo, Fock] = rhf(sys,opts,P_guess)

    Nelec = sys.Nelec;
    Smat = sys.overlap;
    Zmat = sys.e1int;
    VVmat = sys.e2int;
    Vnuc = sys.Vnuc;
    
    diis_size = opts.diis_size;
    level_shift = opts.level_shift;
    maxit = opts.maxit;
    tol = opts.tol;
    
    fprintf('\n==================================++Entering RHF Routine++=================================\n')
    
    Norb = size(Smat,1);    
    Nocc = Nelec/2;

    tic_Start = tic;
    
    % Orthogonalization of basis
    X = orthomat(Smat);
    
    % Initial guess at P 
    if nargin < 3
        P = zeros(Norb);
    else
        P = P_guess;
    end
    
%     % Calculate Fock matrix with guess
    G = einsum_kg(VVmat,P,'prqs,rs->pq') - ...
        0.5*einsum_kg(VVmat,P,'prsq,rs->pq');    
    Fock = Zmat + G;% + level_shift*(Smat-Smat.*P.*Smat);
    
    % Diagonalize Fock and get MO vectors and energies
    [Cp,eps_mo] = eig(X'*Fock*X); eps_mo = diag(eps_mo);
    [~,idx] = sort(eps_mo,'ascend'); eps_mo = eps_mo(idx);   
    C = X*Cp(:,idx);
    
    % Build density matrix
    P = 2*einsum_kg(C(:,1:Nocc),C(:,1:Nocc),'pi,qi->pq');

    Escf = 0.5*einsum_kg(P,Zmat+Fock,'pq,qp->') + Vnuc;

    % DIIS containers
    F_list = zeros(Norb^2,diis_size); F_resid_list = zeros(Norb^2,diis_size);
    
    % SCF Iterations
    it_micro = 0; it_macro = 0; flag_conv = false;
    while it_micro < maxit
        
        tic
        
        % Calculate Fock matrix 
        G = einsum_kg(VVmat,P,'prqs,rs->pq') - ...
            0.5*einsum_kg(VVmat,P,'prsq,rs->pq');

        Fock = Zmat + G + level_shift*(Smat-Smat*P*Smat);

        % Build DIIS Residual
        diis_r = X*(Fock*P*Smat - Smat*P*Fock)*X;
        scf_resid = sqrt(mean(diis_r(:).^2));
      
        % append trial and residual vectors to lists
        F_list(:,mod(it_micro,diis_size)+1) = Fock(:);
        F_resid_list(:,mod(it_micro,diis_size)+1) = diis_r(:);
         
        % diis extrapolate
        if mod(it_micro,diis_size) == 0 && it_micro > 0
           it_macro = it_macro + 1;
           fprintf('\nDIIS Cycle - %d',it_macro)
           Fock = diis_xtrap(F_list,F_resid_list);
           Fock = reshape(Fock,[Norb,Norb]);
        end    
        
        % Diagonalize Fock and get MO vectors and energies
        [Cp,eps_mo] = eig(X'*Fock*X); eps_mo = diag(eps_mo);
        [~,idx] = sort(eps_mo,'ascend'); eps_mo = eps_mo(idx);
        C = X*Cp(:,idx);

        % Compute HF energy and density matrix 
        P = 2*einsum_kg(C(:,1:Nocc),C(:,1:Nocc),'pi,qi->pq');
        Escf = 0.5*einsum_kg(P,Zmat+Fock,'qp,pq->') + Vnuc;

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

% def scf(Nelec,atom_coordinates,Z,Smat,Zmat,VVmat,Vnuc,diis_size_scf,level_shift):
%     
%     print('Beginning SCF routine...')
%     
%     Norb = Smat.shape[0]
%     
%     Nocc0 = int(Nelec/2)
% 
%     t0 = time.time()
%     
%     # Orthogonalization of basis
%     X = orthomat(Smat,0.0,'symmetric')
%     
%     # Initial guess at P 
%     P = np.zeros((Norb,Norb))
%     
%     # Calculate Fock matrix with guess
%     G = np.einsum('prqs,rs->pq',VVmat,P,optimize=True)-\
%         0.5*np.einsum('prsq,rs->pq',VVmat,P,optimize=True)
%     Fock = Zmat + G + level_shift*(Smat-Smat*P*Smat)
%     
%     eps_mo, Cp = np.linalg.eigh(np.dot(X.T,np.dot(Fock,X)))
%     idx = eps_mo.argsort(); eps_mo = eps_mo[idx]
%     C = dot(X,Cp[:,idx])
%     P = 2*np.einsum('pi,qi->pq',C[:,:Nocc0],C[:,:Nocc0],optimize=True)
%     
%     # DIIS containers
%     F_list = []; resid_F = []
%     
%     # SCF Iterations
%     it = 1; flag_conv_scf = False
%     while it < maxit_scf: 
%         
%         # Calculate Fock matrix with guess
%         G = np.einsum('prqs,rs->pq',VVmat,P,optimize=True)-\
%             0.5*np.einsum('prsq,rs->pq',VVmat,P,optimize=True)
%         Fock = Zmat + G + level_shift*(Smat-Smat*P*Smat)
%         
%         # Build DIIS Residual
%         diis_r = X.dot(Fock.dot(P).dot(Smat) - Smat.dot(P).dot(Fock)).dot(X)
%         scf_resid = np.mean(diis_r**2)**0.5
%         
%         if scf_resid <= tol_scf:
%             flag_conv_scf = True
%             it_scf = it
%             break
%             
%         # Append trial & residual vectors to lists
%         # (Note: do NOT use extrapolated matrices in F_list!)
%         if len(F_list) == diis_size_scf:
%             del F_list[0]
%             del resid_F[0]
%         F_list.append(Fock)
%         resid_F.append(diis_r)
%         
%         # DIIS extrapolated Fock matrix
%         if it >= 2:
%             Fock = diis_pulay_solver(F_list, resid_F)
%     
%         # Compute new orbital guess with DIIS Fock matrix
%         eps_mo, Cp = np.linalg.eigh(np.dot(X.T,np.dot(Fock,X)))
%         idx = eps_mo.argsort(); eps_mo = eps_mo[idx]
%         C = dot(X,Cp[:,idx])
%         
%         # Compute HF energy and density matrix 
%         P = 2*np.einsum('pi,qi->pq',C[:,:Nocc0],C[:,:Nocc0],optimize=True)
%         Escf = 0.5*np.einsum('qp,pq->',P,Zmat+Fock,optimize=True) + Vnuc
%         
%         #if flag_scf_verbose:
%         #print('Iter - {}    Residuum = {:>10.12f}    Escf = {:>10.12f}'.format(it, scf_resid, Escf))
%                 
%         it += 1
%         
%     # Compute new orbital guess with DIIS Fock matrix
%     eps_mo, Cp = np.linalg.eigh(np.dot(X.T,np.dot(Fock,X)))
%     idx = eps_mo.argsort(); eps_mo = eps_mo[idx]
%     C = dot(X,Cp[:,idx])
% 
%     # Compute HF energy and density matrix 
%     P = 2*np.einsum('pi,qi->pq',C[:,:Nocc0],C[:,:Nocc0],optimize=True)
% 
%     if flag_conv_scf:
%         print('Completed SCF Routine in {} iterations ({:>.4} seconds)'.format(it_scf, time.time()-t0))
%         print('Erhf = {} Ha\n'.format(Escf))
%         
%     return Escf, C, P, eps_mo, Fock