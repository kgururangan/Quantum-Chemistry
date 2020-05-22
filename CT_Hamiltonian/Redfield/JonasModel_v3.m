close all
clear all
clc

addpath('/Users/karthik/Documents/MATLAB/2DPES/QPT');
addpath('/Users/karthik/Documents/MATLAB/Colormaps');
addpath('/Users/karthik/Documents/MATLAB/Additional_Functions');
addpath('/Users/karthik/Documents/MATLAB/Additional_Functions/regu');

pltc = default_cmap();

n_quanta = 2;
numstates = 3;
zp = 0;

omega_t = 215; % tuning mode
omega_c = 150; % coupling mode

eegap = [800];
bandgap = 12500;
omega_g = [omega_t,omega_c];
n_modes = length(omega_g);
omega_ge = [0, bandgap, bandgap + eegap];
omega_e = omega_g;
N = n_quanta^n_modes;

Dpv = zeros(numstates,n_modes);
% first row of Bv and Dpv MUST be 0's to get the correct ground state block!
Bv(1,:) = ones(1,n_modes);
Dpv(1,:) = zeros(1,n_modes);

drel = linspace(0,1,100);
dcm = 0;
NN = length(drel);
Ncm = length(dcm);
phi = zeros(8,NN,Ncm);
sitepop = zeros(8,NN,Ncm);
Ec = zeros(8,NN,Ncm);

ASEQ = zeros(8,8,NN*Ncm);
count = 1;

for A1 = 1:length(dcm)
    for A2 = 1:length(drel)

Cre = diag(sqrt([1:n_quanta-1]),-1); Ann = Cre';
Num = Cre*Ann;
Dpv = [0, 0;
       dcm(A1), 0;
       drel(A2), 0]; % tuning mode has non-zero displacements, coupling mode does not
J0 = 0;

% ground state
H0_g = zeros(N);
for i = 1:n_modes
    H0temp = omega_g(i)*Num;
    A = 1;
    for j = 1:n_modes
        if i == j
            A = kron(H0temp,A);
        else
            A = kron(eye(n_quanta),A);
        end
    end
    H0_g = H0_g + A;
end

% excited state e1
d1 = Dpv(2,1);
H1_e = H0_g + bandgap*eye(N) - d1*omega_e(1)*kron(eye(n_quanta),Cre+Ann);

% excited state e2
d2 = Dpv(3,1);
H2_e = H0_g + (bandgap+eegap)*eye(N) - d2*omega_e(1)*kron(eye(n_quanta),Cre+Ann);

% linear vibronic coupling term
%H_c = J0*kron(Cre + Ann,eye(n_quanta));
H_c = J0*eye(N);

H_E = kron([1,0;0,0],H1_e)+kron([0,0;0,1],H2_e)+kron([0,1;1,0],H_c);

[V_E,D_E] = eig(H_E);
E_ec = diag(D_E);
[Esort,idx] = sort(E_ec);
V_E = V_E(:,idx);
%
Ne = length(E_ec);

%VV(:,:,
%Ec(:,A2,A1) = E_ec;

ASEQ(:,:,count) = H_E;
count = count+1;


    end
end

% Eigenshuffle
[Vseq,Dseq] = eigenshuffle(ASEQ);
for i = 1:Ne
    Ec(i,:,:) = reshape(Dseq(i,:),NN,Ncm);
    Vj = squeeze(Vseq(:,:,i));
    
end

%% Solve for density matrix evolution
close all
Nt = 300;
c = 3e-2;
dt_ps = 1;
dt = dt_ps*2*pi*c;
t = 0:dt:(Nt-1)*dt;
rho0 = zeros(Ne,Ne);
rho0(1,1) = 1;
L_E = kron(H_E,eye(Ne)) - kron(eye(Ne),H_E);
rhoM = zeros(Ne^2,Nt);
rho = zeros(Ne,Ne,Nt);
for i = 1:Nt
    rhoM(:,i) = expm(-1i*L_E*t(i))*rho0(:);
    rho(:,:,i) = reshape(rhoM(:,i),Ne,Ne);
end
clr = jet(Ne);
hold on;
for i = 1:Ne
    plot(t,squeeze(rho(i,i,:)),'b-.','color',clr(i,:),'Linewidth',2)
    leg{i} = [num2str(i)];
end
hold off;
axis([0,inf,0,1])
legend(leg)
xlabel('t [rad^{-1}cm]')
ylabel('\rho(t)')
grid on
set(gca,'FontSize',17,'Linewidth',2,'Box','off')

%%
figure(2)
for i = 1:4
    subplot(2,2,i)
    contourf(dcm,drel,squeeze(real(Ec(i,:,:))))
    xlabel('d_{CM}')
    ylabel('\Delta d')
    set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
    grid on
    axis square
    colorbar
    %colormap(brewermap([],'PRGn'))
    colormap(cmocean('thermal'))
    title(sprintf('Energy: %u',i))
    
end

figure(3)
for i = 1:4
    subplot(2,2,i)
    contourf(dcm,drel,squeeze(real(phi(i,:,:))))
    xlabel('d_{CM}')
    ylabel('\Delta d')
    set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
    grid on
    axis square
    colorbar
    %colormap(brewermap([],'RdBu'))
    colormap(cmocean('thermal'))
    title(sprintf('phi: %u',i))
end

figure(4)
for i = 1:4
    subplot(2,2,i)
    contourf(dcm,drel,squeeze(real(sitepop(i,:,:))))
    xlabel('d_{CM}')
    ylabel('\Delta d')
    set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
    grid on
    axis square
    colorbar
    %colormap(brewermap([],'RdBu'))
    colormap(cmocean('thermal'))
    title(sprintf('Site: %u',i))
end
 


%% dA = dB

pltc = [pltc; 0 0 0];

[~,II1] = min( abs(drel) );
[~,II2] = min( abs(dcm) );

Ecm = squeeze(Ec(:,II1,:));
Ecm = [fliplr(Ecm),Ecm];
dpl = [fliplr(-dcm),dcm];
phicm = squeeze(sitepop(:,II1,:));

Erel = squeeze(Ec(:,:,II2));
phirel = squeeze(sitepop(:,:,II2));


figure(1)
subplot(221)
hold on;
for i = 1:size(Ec,1)
    h1 = plot(dpl,Ecm(i,:),'b-.','color',pltc(i,:),'Linewidth',2);
    leg{i} = sprintf('%u',i);
end
hold off;
xlabel('d_{CM}')
ylabel('Energy [cm^{-1}]')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
grid on
axis square
ll = legend(leg); set(ll,'FontSize',13,'Location','Best')
%
subplot(222)
hold on;
for i = 1:size(phi,1)
h2 = plot(dcm,phicm(i,:),'k-','color',pltc(i,:),'Linewidth',2);
end
hold off;
ll = legend(leg); set(ll,'FontSize',13,'Location','Best')
xlabel('d_{CM}')
ylabel('\phi [deg]')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
grid on
axis square

subplot(223)
hold on;
for i = 1:size(Ec,1)
    plot(drel,Erel(i,:),'b-.','color',pltc(i,:),'Linewidth',2);
    leg{i} = sprintf('%u',i);
end
hold off;
xlabel('\Delta d')
ylabel('Energy [cm^{-1}]')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
grid on
axis square
ll = legend(leg); set(ll,'FontSize',13,'Location','Best')
%
subplot(224)
hold on;
for i = 1:size(phi,1)
 plot(drel,phirel(i,:),'k-','color',pltc(i,:),'Linewidth',2);
end
hold off;
ll = legend(leg); set(ll,'FontSize',13,'Location','Best')
xlabel('\Delta d')
ylabel('P(-)')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
grid on
axis square
    
   
%% Neglect vibrational momentum -> Treat vibrational coordinates as parameters
% Holstein molecular crystal model
% Two-site electronic level system
close all
clc
clear all

pltc = default_cmap();

sigma_z = [1, 0; 0, -1];
sigma_x = [0, 1; 1, 0];

q = linspace(-5,5,500);
d = 1; % A/m*omega_0^2
J = 0.1;
lambda = -0.2; % lambda1 - lambda2; dJ/dq -> when J = lambda*q, conical intersection appears at point when avoided crossing would have occured
eps1 = 0;
eps2 = 1;

Aseq = zeros(2,2,length(q));

for i = 1:length(q)
    
    Hrel = 1/2*q(i)^2*eye(2) + q(i)*d*sigma_z + [eps2, 0; 0, eps1] - (J-lambda*q(i))*sigma_x;
    
    Aseq(:,:,i) = Hrel;
    
end

[Vseq,Dseq] = eigenshuffle(Aseq);

% minimum separation is indeed 2J
[deltaJ,indd] = min( abs(diff(Dseq)) );
deltaJ
q(indd)


hold on
for i = 1:2
    plot(q,Dseq(i,:),'b-.','color',pltc(i,:),'Linewidth',2);
end
hold off
set(gca,'FontSize',17,'Linewidth',1.5,'Box','off')
xlabel('q_{rel}')
ylabel('E/m\omega_{-}^2')
grid on
legend('|E_2\rangle','|E_1\rangle')





 
