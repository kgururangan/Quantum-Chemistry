clear all
clc
close all

% derivative matching point algorithm works, but maybe we're not picking 
% the best matching point...

verbose = 1;
flag_plot = 1;

Cc = 16.857629206; % hbar^2/2 in am.cm^-1.ang^2

% number of nodes in solution we are searching for
NNc = 3;
E_nodes = zeros(1,NNc); 
for num_nodes = 1:NNc

    mass = 1; % amu
    xe = 1.5; % ang
    beta = 1; % 1/ang
    omega = 1; De = 5000; % cm^-1

    Rmin = 0.9; Rmax = 20;
    Nx = 5000;
    dx = (Rmax-Rmin)/(Nx-1); 
    x = Rmin:dx:(Rmax);

    %U = zeros(1,Nx);
    %U = 1/2*mass*omega^2*(x-xe).^2;
    U = De*(1-exp(-beta*(x-xe))).^2;

    load h2-rhf.mat
    Hatowv =  219474.631370;  
    bohrtoang = 0.529177;
    x = RJ*bohrtoang; U = (Ebind - min(Ebind))*Hatowv; dx = abs(x(2)-x(1)); Nx = length(x);

    k2fcn = @(E) 1/Cc*mass*(E-U);

    %E_ex = (pi^2*mass^2*[1:10].^2)/(2*(Rmax-Rmin)^2);
    %E_ex = omega*([1:10]+1/2);

    % initialize containers
    delta_hist = []; energy_hist = []; 
    dpsiL_hist = []; dpsiR_hist = [];
    psiL = []; psiR = [];

    tol = 1e-8; delta = 10; it = 1; maxit = 2000;

    EL = 0; EH = 20000;
    energy = 0.5*(EL+EH);
    %
    while abs(delta) > tol && it < maxit

        % energy function
        %k2 = k2fcn(energy);

        % nodal check
        [~,num_nodes_psi] = detect_nodes_v2(psiL,'left'); 

        % find solution with desired number of nodes
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
        %[ITP,ntp] = detect_nodes_v2(k2,'right'); M = max(ITP);
        deps = 1e-6*energy; et = [energy-deps, energy, energy+deps]; loss = zeros(1,length(et));
        for j = 1:length(et)
            k2 = k2fcn(et(j)); xL = x(1:M); xR = x(M:end); fn = @(J) 1 + dx^2/12*k2(J);
            psiL = numerov_prop(xL,k2,'left');
            psiR = numerov_prop(xR,k2,'right');
            psiR = psiR./psiR(1); psiL = psiL./psiL(end);       
            %fe(j) = error function of derivative at xm
           %fe(j) = -fn(M+1)*psiR(2) + 10*fn(M)*psiL(end) - fn(M-1)*psiL(end-1) - 12*psiL(end);    
           %fe(j) = (-fn(M+1)*psiR(2) + 10*fn(M)*psiL(end) - fn(M-1)*psiL(end-1)) - 12*psiL(end);    
           loss(j) = (1/dx*(psiL(end)-psiL(end-1)+psiR(1)-psiR(2)));
        end

        % calculate d(loss)/dE
        izv = 1/(2*deps)*(loss(3)-loss(1));
        % newton-raphson update
        delta = -loss(2)/izv;
        energy = energy + delta;

        % store values and update iteration count
        energy_hist(it) = energy; delta_hist(it) = delta;
        dpsiR_hist(it) = 1/dx*(psiR(2)-psiR(1)); dpsiL_hist(it) = 1/dx*(psiL(M)-psiL(M-1));

        if verbose == 1
            fprintf('Iter-{%d}: Nodes = %d, Energy = %4.2f, abs(Delta) = %4.2f\n',...
                        it,num_nodes_psi,energy,abs(delta));
        end
        it = it + 1;

    end


    k2 = k2fcn(energy); xL = x(1:M); xR = x(M:end); 
    psiL = numerov_prop(xL,k2,'left');  psiL_norm = normalize_psi(psiL); dpsiL = diff(psiL)/dx;
    psiR = numerov_prop(xR,k2,'right'); psiR_norm = normalize_psi(psiR); dpsiR = diff(psiR)/dx;
    psiL = psiL./psiL(end); psiR = psiR./psiR(1);

    psi = [psiL(1:end),psiR(2:end)];
    psi = normalize_psi(psi);

    if it < maxit 
        flag_exit = 1;
        fprintf('Numerov successfully converged after %d iterations to within tolerance...\n', it)
        fprintf('Energy = %4.2f, Nodes = %d\n',energy,num_nodes_psi);
    else
        flag_exit = 2;
        fprintf('Warning: energy not converged, max iterations (%d) exceeded...\n', maxit)
        fprintf('Energy = %4.2f, Nodes = %d\n',energy,num_nodes_psi);
    end


    if flag_plot == 1
        % Plotting
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.15, 0.6, 0.6]);

        subplot(221)
        plot(x,U,'k-','Linewidth',2); hold on
        h1 = plot(x,energy*ones(1,Nx),'r-','Linewidth',2);
        ll = legend(h1,sprintf('E = %4.2f cm^{-1}',energy)); set(ll,'FontSize',12,'Location','NorthEast');
        axis tight
        xlabel('x / $\AA$','Interpreter','latex')
        ylabel('$U(x)$ / cm$^{-1}$','Interpreter','latex')
        set(gca,'FontSize',17,'Linewidth',2,'Box','off')

        subplot(222)
        plot(xL,psiL,xR,psiR,'Linewidth',2); hold on;
        scatter([x(M),x(M)],[psiL(M),psiR(1)],60,'MarkerEdgeColor',[0,0.8,0],'MarkerFaceColor',[0,0.8,0]); hold off
        xlabel('x / $\AA$','Interpreter','latex')
        ylabel('$\psi(x)$','Interpreter','latex')
        ll = legend('\psi_L(x)/\psi_L(x_M)','\psi_R(x)/\psi_R(x_M)','Match Point'); set(ll,'FontSize',12,'Location','Best');
        yyaxis right
        plot(x,U,'Linewidth',2)
        ylabel('$U(x)$ / cm$^{-1}$','Interpreter','latex')
        set(gca,'FontSize',17,'Linewidth',2,'Box','off');
        axis tight

        subplot(2,2,3:4)
        plot(x,psi,'Linewidth',2); 
        xlabel('x / $\AA$','Interpreter','latex')
        ylabel('$\psi(x)$','Interpreter','latex')
        yyaxis right
        plot(x,U,'k-','Linewidth',2)
        ylabel('$U(x)$ / cm$^{-1}$','Interpreter','latex')
        set(gca,'FontSize',17,'Linewidth',2,'Box','off');
        axis tight
    end

    E_nodes(num_nodes) = energy;

end