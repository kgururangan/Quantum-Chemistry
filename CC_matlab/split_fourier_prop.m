function [WF,E] = split_fourier_prop(sim_par, Vfcn, flag_imag)

% Split operator fourier transform method
% based on Trotter expansion of unitary time evolution for single particle
% Hamiltonian H = p^2/2m + V(x) 
% U = exp(-iHt) = exp(-i*V*dt/2)exp(-ik''x)exp(-i*k^2/2m*dt)exp(ik'x)exp(-i*V*dt/2)
% U_R = exp(-i*V*dt/2); U_K = exp(-i*k^2/2m*dt)
% <x|psi_t> = U_R*ifft(U_K*fft(U_R*psi_0)))

if nargin < 3
    flag_imag = false;
end

if hasfield(sim_par,'resolution')
    res_t = sim_par.resolution;
else
    res_t = 1;
end

if hasfield(sim_par,'plot')
    flag_plot = sim_par.plot;
else
    flag_plot = false;
end

t = sim_par.t; Nt = length(t); tmax = max(abs(t)); dt = abs(t(2)-t(1));
x = sim_par.x; Nx = length(x); xmax = max(abs(x)); dx = abs(x(2)-x(1));
if sim_par.flag_centered 
    k = [0:1:Nx/2-1,-Nx/2:1:-1]*pi/xmax; 
else
    k = [-Nx/2:1:-1, 0:1:Nx/2-1,]*pi/xmax;
end

H_k = @(k) k.^2/2;
 
if flag_imag
    Ux = @(x) exp(-0.5*Vfcn(x)*dt);
    Uk = @(k) exp(-H_k(k)*dt);
else
    Ux = @(x) exp(-1i*0.5*Vfcn(x)*dt);
    Uk = @(k) exp(-1i*H_k(k)*dt);
end

% initial condition
A = sim_par.init;

% allocate storage arrays
nrec = floor(Nt/res_t);
WF = zeros(Nx,nrec);
E = zeros(1,nrec);

WF(:,1) = A;
E(1) = calculate_energy(A,H_k(k'),Vfcn(x'),dx);
ct = 1;

A = reshape(A,Nx,1);

for i = 1:Nt
    
    % propagate by 1/2 timestep in position
    A = Ux(x').*A;
    % fft to k-space
    A = fft(A,Nx);
    % propagate by full timestep in k-space
    A = Uk(k').*A;
    % ifft back to position
    A = ifft(A,Nx);
    % propagate by 1/2 timestep in position
    A = Ux(x').*A;
    
    % calculate probability density 
    rho = conj(A).*A;
    % normalization factor
    n0 = sum(rho)*dx;
    % normalize WF
    A = A./sqrt(n0);
    
    % store wavefunction
    if mod(i,res_t) == 0
        WF(:,ct) = A;
        E(ct) = calculate_energy(A,H_k(k'),Vfcn(x'),dx);
        ct = ct+1;
    end
    
end

   
if flag_plot
    figure(5)
    for j = 1:nrec
        % plot result
        subplot(211)
        plot(x,sqrt(conj(WF(:,j)).*WF(:,j)),'Linewidth',2)
        xlabel('$x$','Interpreter','latex')
        ylabel('$|\psi(x,t)|$','Interpreter','latex')
        grid on
        set(gca,'FontSize',17,'Linewidth',2,'Box','off')
        ll = legend(sprintf('t = %d dt',j*res_t)); set(ll,'Location','NorthEast');
        axis([-inf,inf,0,max(max(abs(WF)))])
        pause(0.001)

        subplot(212)
        plot(E(1:j),'Linewidth',2)
        xlabel('Timestep/$\Delta t$','Interpreter','latex')
        ylabel('Energy/$\hbar \omega$','Interpreter','latex')
        set(gca,'FontSize',17,'Linewidth',2,'Box','off')
        grid on
        axis([0,nrec,min(E)-1e-6,max(E)+1e-6])
    end
end

end
    


% functions
function Eout = calculate_energy(wfc, H_k, H_r, dx)
    % Creating momentum and conjugate wavefunctions
    wfc_k = fft(wfc);
    wfc_c = conj(wfc);

    % Finding the momentum and real-space energy terms
    energy_k = wfc_c.*ifft(H_k .* wfc_k);
    energy_r = wfc_c.* H_r .* wfc;

    % Integrating over all space
    energy_final = 0;
    for i = 1:length(energy_k)
        energy_final  = energy_final + real(energy_k(i) + energy_r(i));
    end 

    Eout =  energy_final*dx;
end
