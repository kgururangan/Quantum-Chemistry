function [Xfit,Yfit,X0] = holstein_absorption_fit_v2(Ydata,omega,nruns,maxeval,lb,ub,n_quanta,n_modes,nstates,fano,x0)

    Xfit = cell(1,size(Ydata,1));
    Yfit = cell(1,size(Ydata,1));
    X0 = cell(nruns,size(Ydata,1));
    
    options = optimoptions('LSQNONLIN','Display','iter-detailed','Diagnostics','on','Algorithm','trust-region-reflective');
    options.FunctionTolerance = 1e-50;
    options.OptimalityTolerance = 1e-50;
    options.StepTolerance = 1e-50;
    options.FunValCheck = 'On';
    options.MaxFunEvals = maxeval;
    options.MaxIterations = 1e3; 

    for i = 1:size(Ydata,1)

        Xs = cell(1,nruns);
        fv = zeros(1,nruns);

        Y = Ydata(i,:);

        for k = 1:nruns
            fun = @(x) ( ( Y - holstein_abs(x,omega,n_quanta,n_modes,nstates,fano) ) ).^2;
            if isempty(x0)
                x0 = (ub-lb).*rand(1,length(ub))+lb;
            end
            X0{k,i} = x0;
            [xsol,fval,residual,exitflag,~,~,J] = lsqnonlin(fun,x0,lb,ub,options);
            fv(k) = fval;
            Xs{k} = xsol;
        end

        [~,ii] = min(fv);
        xsolve = Xs{ii};

        Xfit{i} = xsolve;
        Yfit{i} = holstein_abs(Xfit{i},omega,n_quanta,n_modes,nstates,fano);
    end
    
end

%
function Yout = holstein_abs(x,omega,n_quanta,n_modes,nstates,fano)

    G = zeros(n_modes,n_modes,n_modes);
    zp = 0;
    kbT = 200;
    pop_min = 0.01;
    pmin = 1e-4;
    N = n_quanta^n_modes;
    
    if fano ~= 1
        % lorentzian, symmetric
        lne = @(w,g,q) 1./(1 + ((omega-w)./g).^2 );
        fano_q = 0;
    else
        % fano, asymmetric
        lne = @(w,g,q) (q +  (omega-w)./g ) .^2./(1 + ((omega-w)./g).^2 );
        fano_q = x(end);
    end
    
    omega_g = x(1:n_modes);
    omega_ge = x(n_modes+1:n_modes+nstates); omega_ge = [0, omega_ge];
    Dpv = x(n_modes+nstates+1:n_modes+nstates+nstates*n_modes);
    Bv = x(n_modes+nstates+nstates*n_modes+1:n_modes+nstates+2*nstates*n_modes);
    J = x(n_modes+nstates+2*nstates*n_modes+1:n_modes+nstates+2*nstates*n_modes+nstates^2);
    gamma = x(n_modes+nstates+2*nstates*n_modes+nstates^2+1);
    
    DPV = reshape(Dpv,nstates,n_modes); DPV = [zeros(1,n_modes);DPV];
    BV = reshape(Bv,nstates,n_modes); BV = [ones(1,n_modes);BV] + (1e-6)*ones(nstates+1,n_modes);
    Jcoul = zeros(nstates+1,nstates+1); Jcoul(2:end,2:end) = reshape(J,nstates,nstates);
    
    [Htot,FCmat] = holstein_hamiltonian(n_quanta,omega_g,omega_ge,nstates+1,BV,DPV,Jcoul,G,zp);
    [V,D] = eig(Htot); E = diag(D); Vt = V';
    Z = trace(expm(-D/kbT));
    rho_eq = expm(-D/kbT)/Z;
    
    Nalph = length(rho_eq(1:N,1:N) >= pop_min);   Ialph = 1:Nalph;
    Ng = N;  Ig = 1:N;
    Ne = length(N+1:length(E));  Ie = N+1:length(E);
    E_alph = E(Ialph); E_g = E(Ig); E_ec = E(Ie);
    
    %V1 = V(Ialph,Ie)
    %V2 = V(Ie,Ialph)
    %Ialph
    %Ie
    %V
    
    % (3a) fermionic operators
    cn = cell(nstates,N);
    e = 1;
    for i = 1:nstates
        g = 1;
        for j = 1:N
            temp = zeros(N*(nstates+1));
            temp(g,e+Ng) = 1;
            cn{i,j} = temp;
            g = g + 1;
            e = e + 1;
        end
    end
    
    % (3b) dipole operator
    % bare electronic dipole moment (ground to each exc manifold)
    mu0 = ones(nstates,1);

    W = 0;
    for i = 1:nstates
        temp = zeros(N*(nstates+1));
        for j = 1:N
            temp = temp + (cn{i,j}+cn{i,j}');
        end
        W = W - mu0(i)*temp;
    end
    
    P2 = zeros(Nalph,Ne);
    for alpha = 1:Nalph
        for a = 1:Ne
            %P2(alpha,a) = rho_eq(alpha,alpha)*V(Ialph(alpha)+N,Ie(a))*Vt(Ialph(alpha)+N,Ie(a));
            P2(alpha,a) = rho_eq(alpha,alpha)*W(Ialph(alpha),Ie(a))*W(Ialph(alpha),Ie(a));
            %P2(alpha,a) = rho_eq(alpha,alpha)*(FCmat(Ialph(alpha),Ie(a))*transpose(FCmat(Ialph(alpha),Ie(a))));
        end
    end
    P2
    ind = find( abs(P2) > pmin); 
    [ialph,ia] = ind2sub(size(P2),ind);
    
    Yout = zeros(1,length(omega));
    for alpha = 1:length(ialph)
        for a = 1:length(ia)
            w0 = E_ec(ia(a))-E_alph(ialph(alpha));
            f2 = P2(ialph(alpha),ia(a))*lne(w0,gamma,fano_q);
            Yout = Yout + f2;
        end
    end 
    
    if isempty(ind)
        Yout = zeros(1,length(omega));
    else
        Yout = Yout/max(Yout);
    end
    
end

