function [energy, psi, energy_hist, delta_hist, flag_exit] = numerov_function_v2(x,U,mass,Elim,num_nodes,maxit,tol,flag_verbose)

    % Units: energy in cm^-1
    %        masses in amu
    %        length in ang

    Cc = 16.857629206; % hbar^2/2 in amu.cm^-1.ang^2
    dx = abs(x(2)-x(1));
    Nx = length(x);

    k2fcn = @(E) 1/Cc*mass*(E-U);

    % initialize containers
    delta_hist = []; energy_hist = []; %dpsiL_hist = []; dpsiR_hist = [];
    psiL = []; psiR = [];

    % convergence and iteration parameters
    %tol = 1e-7; 
    delta = 10; it = 1; 

    % initialize energy bounds for given node count
    if isempty(Elim)
        EL = 0; EH = U(end);
    else
        EL = Elim(1); EH = Elim(2);
    end

    % initial bisection guess at energy
    energy = 0.5*(EL+EH);

    % numerov algorithm loop
    while abs(delta) > tol && it < maxit

        % nodal check
        [~,num_nodes_psi] = detect_nodes_v2(psiL,'left'); 

        % find a solution with desired number of nodes
        while num_nodes_psi ~= num_nodes
            k2 = k2fcn(energy);
            [ITP,ntp] = detect_nodes_v2(k2,'left');
            if ntp == 2 && max(ITP) < Nx
                M = max(ITP);
                xL = x(1:M);
                psiL = numerov_prop(xL,k2,'left');
                [~,num_nodes_psi] = detect_nodes_v2(psiL,'left');
                if num_nodes_psi > num_nodes
                    EH = energy;
                elseif num_nodes_psi < num_nodes
                    EL = energy;
                else
                    break
                end
            else
                EH = 0.99*EH;
                continue
            end
            energy = 0.5*(EH+EL);
        end

        % if nodes are good, match derivatives using Newton Raphson
        deps = 1e-6*energy; et = [energy-deps, energy, energy+deps]; loss = zeros(1,length(et));
        for j = 1:length(et)
            k2 = k2fcn(et(j)); xL = x(1:M); xR = x(M:end); %fn = @(J) 1 + dx^2/12*k2(J);
            psiL = numerov_prop(xL,k2,'left');
            psiR = numerov_prop(xR,k2,'right');
            psiR = psiR./psiR(1); psiL = psiL./psiL(end);       
            loss(j) = (1/dx*(psiL(end)-psiL(end-1)+psiR(1)-psiR(2)));
        end

        % calculate d(loss)/dE
        izv = 1/(2*deps)*(loss(3)-loss(1));

        % newton-raphson update
        delta = -loss(2)/izv;

        if it > 10 && length(unique(delta_hist)) < 4
            energy = energy + alpha_w*delta;
            fprintf('Damping applied...\n')
        else
            energy = energy + delta;
        end

        % store values and update iteration count
        energy_hist(it) = energy; delta_hist(it) = delta;
        %dpsiR_hist(it) = 1/dx*(psiR(2)-psiR(1)); dpsiL_hist(it) = 1/dx*(psiL(M)-psiL(M-1));

        % print to screen if verbose
        if flag_verbose == 1
            fprintf('Iter-{%d}: Nodes = %d, Energy = %4.2f, abs(Delta) = %4.2f\n',...
                        it,num_nodes_psi-1,energy,abs(delta));
        end

        it = it + 1;
    end

    % Print final energy and loop exit condition
    if it < maxit 
        flag_exit = 0;
        fprintf('Numerov successfully converged after %d iterations to within tolerance...\n', it)
        fprintf('Energy = %4.2f, Nodes = %d\n',energy,num_nodes_psi-1);
    else
        flag_exit = 1;
        fprintf('Warning: energy not converged, max iterations (%d) exceeded...\n', maxit)
        fprintf('Energy = %4.2f, Nodes = %d\n',energy,num_nodes_psi-1);
    end

    % Calculate final psi at given node count with final energy
    k2 = k2fcn(energy); xL = x(1:M); xR = x(M:end); 
    psiL = numerov_prop(xL,k2,'left');  %psiL_norm = normalize_psi(psiL); dpsiL = diff(psiL)/dx;
    psiR = numerov_prop(xR,k2,'right'); %psiR_norm = normalize_psi(psiR); dpsiR = diff(psiR)/dx;
    psiL = psiL./psiL(end); psiR = psiR./psiR(1);
    psi = [psiL(1:end),psiR(2:end)];
    psi = normalize_psi(psi);

end

function [psi] = numerov_prop(x,k2,dir)
    %Nx = floor(abs(Rlim(2)-Rlim(1))/dx)+1;
    Nx = length(x); dx = abs(x(2)-x(1)); dx2 = dx^2;
    psi = zeros(1,Nx);
    %dx = abs(Rlim(2)-Rlim(1))/(Nx-1);
    fn = @(i) 1+dx2/12*k2(i); % FIX k2(i) HERE!!! NOT THE RIGHT ONE FOR LEFT AND RIGHT BOUNDARIES!!!
    if strcmp(dir,'left')
        psi(1) = 0; psi(2) = 1e-20;
        for i = 2:Nx-1
            psi(i+1) = ( (12-10*fn(i))*psi(i) - fn(i-1)*psi(i-1) )/fn(i+1);
        end
    else
        psi(end) = 1e-30; % arbitrary small number > 0 since we are not at +infinity
        psi(end-1) = exp(x(end-1)*sqrt(-k2(end-1)))/exp(x(end)*sqrt(-k2(end)))*psi(end); % JWKB formula
        for i = Nx-1:-1:2
            psi(i-1) = ( (12-10*fn(i))*psi(i) - fn(i+1)*psi(i+1) )/fn(i-1);
        end

    end
end

function [ZC,NN] = detect_nodes_v2(y,dir)    
    if isempty(y)
        ZC = []; NN = 1e3;
        return 
    end
    NN = 0; ZC = [];
    if strcmp(dir,'left')
        for i = 1:length(y)-1
            if (y(i) > 0 && y(i+1) <0) || (y(i) < 0 && y(i+1) > 0) || (y(i) == 0)
                NN = NN + 1;
                ZC(NN) = i;
            end
        end
    else
        for i = 2:length(y)
            if (y(i) > 0 && y(i-1) <0) || (y(i) < 0 && y(i-1) > 0) || (y(i) == 0)
                NN = NN + 1;
                ZC(NN) = i;
            end
        end
    end   
end

function [ psi_norm ] = normalize_psi( psi )
    psi_norm = psi./(sqrt(sum(conj(psi).*psi)));
end