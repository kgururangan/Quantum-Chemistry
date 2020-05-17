function [t1,t2,Ecorr] = ccsd(sys,opts)

    VM = sys.VM; FM = sys.FM; occ = sys.occ; unocc = sys.unocc;

    diis_size = opts.diis_size;
    maxit = opts.maxit;
    tol = opts.tol;
    
    tic
    fprintf('\n==============++Entering CCSD Routine++=================\n')
    
    Nocc = length(occ);
    Nunocc = length(unocc);
    Nov = Nocc*Nunocc; Noovv = Nov^2; Hbar_dim = Nov + Noovv;
    
    % Initialize T1 and T2 DIIS containers
    T = zeros(Hbar_dim,1);
    T_list = zeros(Hbar_dim,diis_size);
    T_resid_list = zeros(Hbar_dim,diis_size);
    
    % Jacobi/DIIS iterations
    it = 1; flag_conv = 0;
    while it < maxit
        
        % store old T and get current diis dimensions
        T_old = T;
        
        % build t1 and t2 diagrammatically
        [X_ai, X_abij] = build_t1_t2(reshape(T(1:Nov),Nunocc,Nocc),reshape(T(Nov+1:end),Nunocc,Nunocc,Nocc,Nocc),...
                                    VM,FM,occ,unocc);
                                
        % update t1 and t2 by Jacobi                        
        [t1,t2] = update_t1_t2(X_ai,X_abij,VM,FM,occ,unocc);
        
        % store vectorized results
        T(1:Nov) = t1(:); T(Nov+1:end) = t2(:);

        % build DIIS residual
        T_resid = T - T_old;
        
        % check for exit condition
        ccsd_resid = sqrt(mean(T_resid.^2));
        if ccsd_resid < tol
            flag_conv = 1;
            break;
        end
        
        % append trial and residual vectors to lists
        T_list(:,mod(it,diis_size)+1) = T;
        T_resid_list(:,mod(it,diis_size)+1) = T_resid;
         
        % diis extrapolate
        if it >= 2*diis_size
           T = diis_xtrap(T_list,T_resid_list);
        end
        
        t1 = reshape(T(1:Nov),Nunocc,Nocc);
        t2 = reshape(T(Nov+1:end),Nunocc,Nunocc,Nocc,Nocc);

        % CC correlation energy
        Ecorr = cc_energy(t1,t2,VM,FM,occ,unocc);
        
        fprintf('\nIter-%d     Residuum = %4.12f      Ecorr = %4.12f',it,ccsd_resid,Ecorr);
        
        it = it + 1;
        
    end
    
    
    if flag_conv == 1
        fprintf('\nCCSD successfully converged in %d iterations (%4.2f seconds)\n',it,toc);
    else
        fprintf('\nCCSD failed to converged in %d iterations\n',maxit)
    end
            
        

end


       