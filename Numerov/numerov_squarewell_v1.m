clear all
clc
close all

Cc = 16.857629206; % hbar^2/2 in am.cm^-1.ang^2

Rmin = 0.9; Rmax = 20;
Nx = 5000;
dx = (Rmax-Rmin)/(Nx-1); 
mass = 1;

xe = 2; beta = 1; omega = 1; De = 20;
x = Rmin:dx:(Rmax);
%U = zeros(1,Nx);
%U = 1/2*mass*omega^2*(x-xe).^2;
U = De*(1-exp(-beta*(x-xe))).^2;
M = floor(Nx/2);

%E_ex = (pi^2*mass^2*[1:10].^2)/(2*(Rmax-Rmin)^2);
E_ex = omega*([1:10]+1/2);

k2fcn = @(E) 2*mass*(E-U);

EL = 0; EH = 10;

tol = 1e-7; delta = 10; it = 1; maxit = 500;
delta_hist = []; energy_hist = []; dpsiL_hist = []; dpsiR_hist = [];
psiL = []; psiR = [];

num_nodes = 2;
energy = 0.5*(EL+EH);
while abs(delta) > tol && it < maxit
    

    k2 = k2fcn(energy);
    
    [~,num_nodes_psi] = detect_nodes_v2(psiL,'left'); 
    
    % nodal detection
    while num_nodes_psi ~= num_nodes
        k2 = k2fcn(energy);
        [ITP,ntp] = detect_nodes_v2(k2,'right');
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
    deps = 1e-4*energy; et = [energy-deps, energy, energy+deps]; fe = zeros(1,length(et));
    for j = 1:length(et)
        k2 = k2fcn(et(j)); xL = x(1:M); xR = x(M:end); fn = @(J) 1 + dx^2/12*k2(J);
        psiL = numerov_prop(xL,k2,'left');
        psiR = numerov_prop(xR,k2,'right');
        psiR = psiR./psiR(1); psiL = psiL./psiL(end);       
        %fe(j) = error function of derivative at xm
       %fe(j) = -fn(M+1)*psiR(2) + 10*fn(M)*psiL(end) - fn(M-1)*psiL(end-1) - 12*psiL(end);    
       %fe(j) = (-fn(M+1)*psiR(2) + 10*fn(M)*psiL(end) - fn(M-1)*psiL(end-1)) - 12*psiL(end);    
       fe(j) = (1/dx*(psiL(end)-psiL(end-1)+psiR(1)-psiR(2)))^2;
    end
    dpsiR_hist(it) = 1/dx*(psiR(2)-psiR(1)); dpsiL_hist(it) = 1/dx*(psiL(M)-psiL(M-1));
    izv = 1/(2*deps)*(fe(3)-fe(1));
    delta = -fe(2)/izv;
    energy = energy + delta;
    
    energy_hist(it) = energy; delta_hist(it) = delta;
    
    fprintf('Iter-{%d}: Nodes = %d, Energy = %4.2f, abs(Delta) = %4.2f\n',...
                it,num_nodes_psi,energy,abs(delta));
    
    it = it + 1;
    
end

% Plotting
psi = [psiL(1:end),psiR(2:end)];
psi = normalize_psi(psi);

subplot(221)
plot(x,U,x,energy*ones(1,Nx),'k-','Linewidth',2)
axis tight

subplot(222)
plot(xL,normalize_psi(psiL),xR,normalize_psi(psiR))
yyaxis right
plot(x,U,'Linewidth',2)
axis tight

subplot(2,2,3:4)
plot(x,psi,'Linewidth',2);
yyaxis right
plot(x,U,'Linewidth',2)
axis tight
            
            %%

    
    
    
    
    
    
    
    
    
    
    
    
    
           % setting them equal to 1 enforces even nodes
        psiL = psiL/psiL(M);
        psiR = psiR/psiR(1);
        % setting them to 0 should enforce odd but how
        %psiL(M) = 0; psiR(1) = 0;
        dpsiL = 1/dx*(psiL(M)-psiL(M-1));
        dpsiR = 1/dx*(psiR(2)-psiR(1));
        delta = dpsiR - dpsiL;

    
    delta_hist(it) = delta; energy_hist(it) = energy;
    [~,nodesL] = detect_nodes_v2(psiL,'left'); [~,nodesR] = detect_nodes_v2(psiR,'right');
    nodes_psi = nodesL+nodesR-2;
    range_E = EH-EL;
    fprintf('Iter-{%d}: Nodes = %d, Energy = %4.2f, abs(Delta) = %4.2f, Erange = %4.2f\n',it,nodes_psi,energy, abs(delta), range_E);
    it = it + 1;
    
    
%     if delta-delta0 < 1e-9
%         EH = 1.1*EH;
%     end
%     
%end

psi = [psiL(1:end),psiR(2:end)];
psi = normalize_psi(psi);

subplot(221)
plot(x,U,x,energy*ones(1,Nx),'k-','Linewidth',2)
axis tight

subplot(222)
plot(xL,normalize_psi(psiL),xR,normalize_psi(psiR))
yyaxis right
plot(x,U,'Linewidth',2)
axis tight

subplot(2,2,3:4)
plot(x,psi,'Linewidth',2);
yyaxis right
plot(x,U,'Linewidth',2)
axis tight

%%
clear all
clc
close all

Cc = 16.857629206;

dx = 0.001;
Nx = 1024;
x = 0:dx:(Nx-1)*dx;
U = 1000*ones(1,Nx);

EL = 0; EH = 500;
energy = 0.5*(EH+EL);
k2 = 1/Cc*(energy-U);

psi = numerov_prop(x,k2,'left');
plot(x,psi)

