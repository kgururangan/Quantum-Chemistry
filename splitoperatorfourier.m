clear all
clc
close all

% Split operator fourier transform method
% based on Trotter expansion of unitary time evolution for single particle
% Hamiltonian H = p^2/2m + V(x) 
% U = exp(-iHt) = exp(-i*V*dt/2)exp(-ik''x)exp(-i*k^2/2m*dt)exp(ik'x)exp(-i*V*dt/2)
% U_R = exp(-i*V*dt/2); U_K = exp(-i*k^2/2m*dt)
% <x|psi_t> = U_R*ifft(U_K*fft(U_R*psi_0)))

imag_time = 0;

Nt = 128;
Nx = 1024; 
xmax = 10;
tmax = 5;

x = linspace(-xmax,xmax,Nx); dx = abs(x(2)-x(1));
% the key is in setting the frequency k-space axis corretly!
%k = [0:1:Nx/2-1,-Nx/2:1:-1]*pi/xmax; %dk = abs(k(2)-k(1));
%k = ifftshift([-Nx/2:1:Nx/2-1]*pi/xmax);
k = ifftshift(FFTaxis_kg(x,1,0));
t = linspace(0,tmax,Nt); dt = abs(t(2)-t(1));
%
H_k = @(k) k.^2/2;
 
% harmonic trap
x0 = 0.5;
k0 = 2;
V = @(x) 1/2*k0*(x-x0).^2;

% notch potential
% x0 = 0.5;
% k0 = 10;
% V = @(x) k0*heaviside(x-x0) + k0*heaviside(x0-4-x);

% delta function potential
% x0 = 0.5;
% [~,i0] = min(abs(x-x0));
% k0 = -1;
% V = @(x) k0*deltafcn(x,x(i0));

% coulomb potential
% k0 = 1;
% V = @(x) k0*exp(-abs(x))./abs(x);


if imag_time == 1
    Ux = @(x) exp(-0.5*V(x)*dt);
    Uk = @(k) exp(-H_k(k)*dt);
else
    Ux = @(x) exp(-1i*0.5*V(x)*dt);
    Uk = @(k) exp(-1i*H_k(k)*dt);
end

% initial condition
sigma0 = 1;
A = 1/sqrt(2*pi*sigma0^2)*exp(-x.^2/(2*sigma0^2));
%A = rand(1,Nx);
%A = ones(1,Nx);
%A = cos(2*pi*x);/

% allocate storage arrays
res_t = 1;
nrec = floor(Nt/res_t);
WF = zeros(Nx,nrec);
E = zeros(1,nrec);

WF(:,1) = A;
E(1) = calculate_energy(A,H_k(k'),V(x'),dx);
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
        E(ct) = calculate_energy(A,H_k(k'),V(x'),dx);
        ct = ct+1;
    end
    
end

%% nonlinear potentials

imag_time = 0;

Nt = 128;
Nx = 1024; 
xmax = 10;
tmax = 50;

x = linspace(-xmax,xmax,Nx); dx = abs(x(2)-x(1));
% the key is in setting the frequency k-space axis corretly!
%k = [0:1:Nx/2-1,-Nx/2:1:-1]*pi/xmax; dk = abs(k(2)-k(1));
k = ifftshift([-Nx/2:1:Nx/2-1]*pi/xmax);
t = linspace(0,tmax,Nt); dt = abs(t(2)-t(1));

H_k = @(k) k.^2/2;

% initial condition
A = exp(-x.^2/2);
%A = cos(2*pi*x);
%A = exp(1i*x).*sech(x);

% allocate storage arrays
res_t = 1;
nrec = floor(Nt/res_t);
WF = zeros(Nx,nrec);
E = zeros(1,nrec);

WF(:,1) = A;
E(1) = calculate_energy(A,H_k(k'),abs(A).*abs(A),dx);
ct = 1;

A = reshape(A,Nx,1);

k1 = 1; k2 = 10;
for i = 1:Nt
    
    V = @(x) k1*abs(A).^2 + k2*abs(A).^4;
    
    if imag_time == 1
        Ux = @(x) exp(-0.5*V(x)*dt);
        Uk = @(k) exp(-H_k(k)*dt);
    else
        Ux = @(x) exp(-1i*0.5*V(x)*dt);
        Uk = @(k) exp(-1i*H_k(k)*dt);
    end
    
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
        E(ct) = calculate_energy(A,H_k(k'),V(x'),dx);
        ct = ct+1;
    end
    
end

    
%%
close all

figure(1)
for j = 1:nrec
    % plot result
    subplot(211)
    plot(x,sqrt(conj(WF(:,j)).*WF(:,j)),x,V(x),'Linewidth',2)
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
    
%%

clear all
clc
close all

dt = 0.001;
fs = 1/dt;
Nt = 1024;
t = [0:(Nt-1)]*dt;

f0 = [5,22,31]; % Hz
f = zeros(1,Nt);
for i = 1:length(f0)
    f = f + cos(2*pi*f0(i)*t);% + rand(1,Nt);
end

%F = fftshift(fft(f,Nt));
F = fft(f,Nt);

% for symmetric fft (real signals)
%F = F(1:floor(Nt/2));

omega = FFTaxis_kg(t,1,1,0);

plot(omega,fftshift(abs(F)))


%% functions

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

function [y] = deltafcn(x,x0)
    y = zeros(length(x),1);
    [~,i0] = min(abs(x-x0));
    y(i0) = 1;
%     if x == x0
%         y = 1;
%     else
%         y = 0;
%     end
end
    
    