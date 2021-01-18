function [cc_t,Ecorr] = uccsdt2_full(sys,opts,T_guess)

    diis_size = opts.diis_size;
    maxit = opts.maxit;
    tol = opts.tol;
    shift = opts.shift;

    if nargin < 3 || isempty(T_guess)
        T = zeros(sys.triples_dim,1);
    else
        T = T_guess;
    end
    
    tic_Start = tic;
    fprintf('\n==================================++Entering UCCSDT-2 Routine++=================================\n')

    % Initialize T1, T2, and T3 DIIS containers
    
    T_list = zeros(sys.triples_dim,diis_size);
    T_resid_list = zeros(sys.triples_dim,diis_size);

    szt1a = [sys.Nvir_alpha, sys.Nocc_alpha];
    szt1b = [sys.Nvir_beta, sys.Nocc_beta];
    szt2a = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nocc_alpha, sys.Nocc_alpha];
    szt2b = [sys.Nvir_alpha, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_beta];
    szt2c = [sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_beta, sys.Nocc_beta];
    szt3a = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nvir_alpha, sys.Nocc_alpha, sys.Nocc_alpha, sys.Nocc_alpha];
    szt3b = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_alpha, sys.Nocc_beta];
    szt3c = [sys.Nvir_alpha, sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_beta, sys.Nocc_beta];
    szt3d = [sys.Nvir_beta, sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_beta, sys.Nocc_beta, sys.Nocc_beta];

    % Jacobi/DIIS iterations
    it_micro = 0; flag_conv = 0; it_macro = 0; Ecorr_old = 0.0;
    %fprintf('\nDIIS Cycle - %d',it_macro)
    while it_micro < maxit
        
        tic
        
        % store old T and get current diis dimensions
        T_old = T;
        
        t1a = reshape(T(sys.posv{1}),szt1a);
        t1b = reshape(T(sys.posv{2}),szt1b);
        t2a = reshape(T(sys.posv{3}),szt2a);
        t2b = reshape(T(sys.posv{4}),szt2b);
        t2c = reshape(T(sys.posv{5}),szt2c);
        t3a = reshape(T(sys.posv{6}),szt3a);
        t3b = reshape(T(sys.posv{7}),szt3b);
        t3c = reshape(T(sys.posv{8}),szt3c);
        t3d = reshape(T(sys.posv{9}),szt3d);

        % CC correlation energy
        Ecorr = ucc_energy(t1a,t1b,t2a,t2b,t2c,sys);
       
        % update t1 and t2 by Jacobi                          
        t1a = update_t1a_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        t1b = update_t1b_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        t2a = update_t2a_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        t2b = update_t2b_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        t2c = update_t2c_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);

        t3a = update_t3a(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        [t3a] = zero_t3a_outside(t3a,sys);

        t3b = update_t3b(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        [t3b] = zero_t3b_outside(t3b,sys);

        t3c = update_t3c(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        [t3c] = zero_t3c_outside(t3c,sys);

        t3d = update_t3d(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        [t3d] = zero_t3d_outside(t3d,sys);


        % store vectorized results
        T(sys.posv{1}) = t1a(:); T(sys.posv{2}) = t1b(:);
        T(sys.posv{3}) = t2a(:); T(sys.posv{4}) = t2b(:); T(sys.posv{5}) = t2c(:);
        T(sys.posv{6}) = t3a(:); T(sys.posv{7}) = t3b(:); T(sys.posv{8}) = t3c(:); T(sys.posv{9}) = t3d(:);

        % build DIIS residual
        T_resid = T - T_old;  
        
        % change in Ecorr
        deltaE = Ecorr - Ecorr_old;

        % check for exit condition
        ccsdt_resid = sqrt(mean(T_resid.^2));
        if ccsdt_resid < tol && abs(deltaE) < tol
            flag_conv = 1;
            break;
        end

        % append trial and residual vectors to lists
        T_list(:,mod(it_micro,diis_size)+1) = T;
        T_resid_list(:,mod(it_micro,diis_size)+1) = T_resid;
         
        % diis extrapolate
        if mod(it_micro,diis_size) == 0 && it_micro > 1
           it_macro = it_macro + 1;
           fprintf('\nDIIS Cycle - %d',it_macro)
           T = diis_xtrap(T_list,T_resid_list);
        end        

        % extract time per iteration in minutes and seconds
        toc_S = toc; toc_M = floor(toc_S/60); toc_S = toc_S - toc_M*60;

        fprintf('\n   Iter-%d    Residuum = %4.10f    dE = %4.10f    Ecorr = %4.10f   (%dm %.1fs)',it_micro,ccsdt_resid,deltaE,Ecorr,toc_M,toc_S);
        
        it_micro = it_micro + 1;
        Ecorr_old = Ecorr;
         
    end

    % return final amplitudes and correlation energy
    cc_t.t1a = reshape(T(sys.posv{1}),szt1a);
    cc_t.t1b = reshape(T(sys.posv{2}),szt1b);
    cc_t.t2a = reshape(T(sys.posv{3}),szt2a);
    cc_t.t2b = reshape(T(sys.posv{4}),szt2b);
    cc_t.t2c = reshape(T(sys.posv{5}),szt2c);
    cc_t.t3a = reshape(T(sys.posv{6}),szt3a);
    cc_t.t3b = reshape(T(sys.posv{7}),szt3b);
    cc_t.t3c = reshape(T(sys.posv{8}),szt3c);
    cc_t.t3d = reshape(T(sys.posv{9}),szt3d);
        
    Ecorr = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);
    
    if flag_conv == 1
        fprintf('\nUCCSDT-2 successfully converged in %d iterations (%4.2f seconds)\n',it_micro,toc(tic_Start));
        fprintf('Total Energy = %4.12f Eh     Ecorr = %4.12f Eh\n',Ecorr+sys.Escf,Ecorr);
    else
        fprintf('\nUCCSDT-2 failed to converged in %d iterations\n',maxit)
    end   
        

end

function [t3a] = zero_t3a_outside(t3a,sys)

    ctP = 0;
    ctQ = 0;

    Nunocc_a = sys.Nvir_alpha;
    Nocc_a = sys.Nocc_alpha;
    
    act_h_alpha_range = [sys.Nocc_alpha - sys.Nact_h_alpha + 1 : sys.Nocc_alpha];
    act_p_alpha_range = [1 : sys.Nact_p_alpha];
    
    for a = 1:Nunocc_a
        for b = a+1:Nunocc_a
            for c = b+1:Nunocc_a
                for i = 1:Nocc_a
                    for j = i+1:Nocc_a
                        for k = j+1:Nocc_a
                            num_act_h = num_active([i,j,k],act_h_alpha_range);
                            num_act_p = num_active([a,b,c],act_p_alpha_range);
                            if num_act_h >=2 && num_act_p >= 2
                                ctP = ctP + 1;
                                continue
                            else          
                                t3a(a,b,c,i,j,k) = 0;
                                t3a(a,b,c,k,i,j) = t3a(a,b,c,i,j,k);
                                t3a(a,b,c,j,k,i) = t3a(a,b,c,i,j,k);
                                t3a(a,b,c,i,k,j) = -t3a(a,b,c,i,j,k);
                                t3a(a,b,c,j,i,k) = -t3a(a,b,c,i,j,k);
                                t3a(a,b,c,k,j,i) = -t3a(a,b,c,i,j,k);

                                % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3a(b,a,c,i,j,k) = -t3a(a,b,c,i,j,k);
                                t3a(b,a,c,k,i,j) = -t3a(a,b,c,i,j,k);
                                t3a(b,a,c,j,k,i) = -t3a(a,b,c,i,j,k);
                                t3a(b,a,c,i,k,j) = t3a(a,b,c,i,j,k);
                                t3a(b,a,c,j,i,k) = t3a(a,b,c,i,j,k);
                                t3a(b,a,c,k,j,i) = t3a(a,b,c,i,j,k);

                                % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3a(a,c,b,i,j,k) = -t3a(a,b,c,i,j,k);
                                t3a(a,c,b,k,i,j) = -t3a(a,b,c,i,j,k);
                                t3a(a,c,b,j,k,i) = -t3a(a,b,c,i,j,k);
                                t3a(a,c,b,i,k,j) = t3a(a,b,c,i,j,k);
                                t3a(a,c,b,j,i,k) = t3a(a,b,c,i,j,k);
                                t3a(a,c,b,k,j,i) = t3a(a,b,c,i,j,k);

                                % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3a(c,b,a,i,j,k) = -t3a(a,b,c,i,j,k);
                                t3a(c,b,a,k,i,j) = -t3a(a,b,c,i,j,k);
                                t3a(c,b,a,j,k,i) = -t3a(a,b,c,i,j,k);
                                t3a(c,b,a,i,k,j) = t3a(a,b,c,i,j,k);
                                t3a(c,b,a,j,i,k) = t3a(a,b,c,i,j,k);
                                t3a(c,b,a,k,j,i) = t3a(a,b,c,i,j,k);

                                % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3a(b,c,a,i,j,k) = t3a(a,b,c,i,j,k);
                                t3a(b,c,a,k,i,j) = t3a(a,b,c,i,j,k);
                                t3a(b,c,a,j,k,i) = t3a(a,b,c,i,j,k);
                                t3a(b,c,a,i,k,j) = -t3a(a,b,c,i,j,k);
                                t3a(b,c,a,j,i,k) = -t3a(a,b,c,i,j,k);
                                t3a(b,c,a,k,j,i) = -t3a(a,b,c,i,j,k);

                                % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3a(c,a,b,i,j,k) = t3a(a,b,c,i,j,k);
                                t3a(c,a,b,k,i,j) = t3a(a,b,c,i,j,k);
                                t3a(c,a,b,j,k,i) = t3a(a,b,c,i,j,k);
                                t3a(c,a,b,i,k,j) = -t3a(a,b,c,i,j,k);
                                t3a(c,a,b,j,i,k) = -t3a(a,b,c,i,j,k);
                                t3a(c,a,b,k,j,i) = -t3a(a,b,c,i,j,k);
                                ctQ = ctQ + 1;
                            end
                            
                        end
                    end
                end
            end
        end
    end
end

function [t3b] = zero_t3b_outside(t3b,sys)

    ctP = 0;
    ctQ = 0;

    Nunocc_a = sys.Nvir_alpha;
    Nunocc_b = sys.Nvir_beta;
    Nocc_a = sys.Nocc_alpha;
    Nocc_b = sys.Nocc_beta;
    
    act_h_alpha_range = [sys.Nocc_alpha - sys.Nact_h_alpha + 1 : sys.Nocc_alpha];
    act_h_beta_range = [sys.Nocc_beta - sys.Nact_h_beta + 1 : sys.Nocc_beta];
    act_p_alpha_range = [1 : sys.Nact_p_alpha];
    act_p_beta_range = [1 : sys.Nact_p_beta];
    
    for a = 1:Nunocc_a
        for b = a+1:Nunocc_a
            for c = 1:Nunocc_b
                for i = 1:Nocc_a
                    for j = i+1:Nocc_a
                        for k = 1:Nocc_b
                            num_act_h = num_active([i,j],act_h_alpha_range) + num_active([k],act_h_beta_range);
                            num_act_p = num_active([a,b],act_p_alpha_range) + num_active([c],act_p_beta_range);

                            if num_act_h >=2 && num_act_p >= 2
                                ctP = ctP + 1;
                                continue
                            else          
                                t3b(a,b,c,i,j,k) = 0;
                                t3b(b,a,c,i,j,k) = -t3b(a,b,c,i,j,k);
                                t3b(a,b,c,j,i,k) = -t3b(a,b,c,i,j,k);
                                t3b(b,a,c,j,i,k) = t3b(a,b,c,i,j,k);
                                ctQ = ctQ + 1;
                            end
                            
                        end
                    end
                end
            end
        end
    end
end

function [t3c] = zero_t3c_outside(t3c,sys)

    ctP = 0;
    ctQ = 0;

    Nunocc_a = sys.Nvir_alpha;
    Nunocc_b = sys.Nvir_beta;
    Nocc_a = sys.Nocc_alpha;
    Nocc_b = sys.Nocc_beta;
    
    act_h_alpha_range = [sys.Nocc_alpha - sys.Nact_h_alpha + 1 : sys.Nocc_alpha];
    act_h_beta_range = [sys.Nocc_beta - sys.Nact_h_beta + 1 : sys.Nocc_beta];
    act_p_alpha_range = [1 : sys.Nact_p_alpha];
    act_p_beta_range = [1 : sys.Nact_p_beta];
    
    for a = 1:Nunocc_a
        for b = 1:Nunocc_b
            for c = b+1:Nunocc_b
                for i = 1:Nocc_a
                    for j = 1:Nocc_b
                        for k = j+1:Nocc_b

                            num_act_h = num_active([i],act_h_alpha_range) + num_active([j,k],act_h_beta_range);
                            num_act_p = num_active([a],act_p_alpha_range) + num_active([b,c],act_p_beta_range);

                            if num_act_h >=2 && num_act_p >= 2
                                ctP = ctP + 1;
                                continue
                            else          
                                t3c(a,b,c,i,j,k) = 0;
                                t3c(a,c,b,i,j,k) = -t3c(a,b,c,i,j,k);
                                t3c(a,b,c,i,k,j) = -t3c(a,b,c,i,j,k);
                                t3c(a,c,b,i,k,j) = t3c(a,b,c,i,j,k);
                                ctQ = ctQ + 1;
                            end
                            
                        end
                    end
                end
            end
        end
    end
end

function [t3d] = zero_t3d_outside(t3d,sys)

    ctP = 0;
    ctQ = 0;

    Nunocc_b = sys.Nvir_beta;
    Nocc_b = sys.Nocc_beta;
    
    act_h_beta_range = [sys.Nocc_beta - sys.Nact_h_beta + 1 : sys.Nocc_beta];
    act_p_beta_range = [1 : sys.Nact_p_beta];
    
    
    for a = 1:Nunocc_b
        for b = a+1:Nunocc_b
            for c = b+1:Nunocc_b
                for i = 1:Nocc_b
                    for j = i+1:Nocc_b
                        for k = j+1:Nocc_b

                            num_act_h = num_active([i,j,k],act_h_beta_range);
                            num_act_p = num_active([a,b,c],act_p_beta_range);

                            if num_act_h >=2 && num_act_p >= 2
                                ctP = ctP + 1;
                                continue
                            else          
                                t3d(a,b,c,i,j,k) = 0;
                                t3d(a,b,c,k,i,j) = t3d(a,b,c,i,j,k);
                                t3d(a,b,c,j,k,i) = t3d(a,b,c,i,j,k);
                                t3d(a,b,c,i,k,j) = -t3d(a,b,c,i,j,k);
                                t3d(a,b,c,j,i,k) = -t3d(a,b,c,i,j,k);
                                t3d(a,b,c,k,j,i) = -t3d(a,b,c,i,j,k);

                                % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3d(b,a,c,i,j,k) = -t3d(a,b,c,i,j,k);
                                t3d(b,a,c,k,i,j) = -t3d(a,b,c,i,j,k);
                                t3d(b,a,c,j,k,i) = -t3d(a,b,c,i,j,k);
                                t3d(b,a,c,i,k,j) = t3d(a,b,c,i,j,k);
                                t3d(b,a,c,j,i,k) = t3d(a,b,c,i,j,k);
                                t3d(b,a,c,k,j,i) = t3d(a,b,c,i,j,k);

                                % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3d(a,c,b,i,j,k) = -t3d(a,b,c,i,j,k);
                                t3d(a,c,b,k,i,j) = -t3d(a,b,c,i,j,k);
                                t3d(a,c,b,j,k,i) = -t3d(a,b,c,i,j,k);
                                t3d(a,c,b,i,k,j) = t3d(a,b,c,i,j,k);
                                t3d(a,c,b,j,i,k) = t3d(a,b,c,i,j,k);
                                t3d(a,c,b,k,j,i) = t3d(a,b,c,i,j,k);

                                % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3d(c,b,a,i,j,k) = -t3d(a,b,c,i,j,k);
                                t3d(c,b,a,k,i,j) = -t3d(a,b,c,i,j,k);
                                t3d(c,b,a,j,k,i) = -t3d(a,b,c,i,j,k);
                                t3d(c,b,a,i,k,j) = t3d(a,b,c,i,j,k);
                                t3d(c,b,a,j,i,k) = t3d(a,b,c,i,j,k);
                                t3d(c,b,a,k,j,i) = t3d(a,b,c,i,j,k);

                                % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3d(b,c,a,i,j,k) = t3d(a,b,c,i,j,k);
                                t3d(b,c,a,k,i,j) = t3d(a,b,c,i,j,k);
                                t3d(b,c,a,j,k,i) = t3d(a,b,c,i,j,k);
                                t3d(b,c,a,i,k,j) = -t3d(a,b,c,i,j,k);
                                t3d(b,c,a,j,i,k) = -t3d(a,b,c,i,j,k);
                                t3d(b,c,a,k,j,i) = -t3d(a,b,c,i,j,k);

                                % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                t3d(c,a,b,i,j,k) = t3d(a,b,c,i,j,k);
                                t3d(c,a,b,k,i,j) = t3d(a,b,c,i,j,k);
                                t3d(c,a,b,j,k,i) = t3d(a,b,c,i,j,k);
                                t3d(c,a,b,i,k,j) = -t3d(a,b,c,i,j,k);
                                t3d(c,a,b,j,i,k) = -t3d(a,b,c,i,j,k);
                                t3d(c,a,b,k,j,i) = -t3d(a,b,c,i,j,k);
                                ctQ = ctQ + 1;
                            end
                            
                        end
                    end
                end
            end
        end
    end
end

function [val] = num_active(x,act_range)
    val = length(find(x >= act_range(1) & x <= act_range(end)));
end

function [flag] = is_P(x,act_range,nact_min)

    nact = length(find(x >= act_range(1) & x <= act_range(end)));
    if nact >= nact_min
        flag = true;
    else
        flag = false;
    end

end