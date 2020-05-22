clear all
close all
clc

% the driving principle behind SF in this model is a lambda or a sequential
% configuration of energies

num_chromophore = 2; % can (feasibly..?) do up to 4 chromophores without crashing
chromophore_spacing = 1;
omega_ge = 12500; % cm^-1
effective_chromophore_charge = ones(1,num_chromophore);

spin_cons = 1;
dets = get_ci_vecs(num_chromophore,spin_cons);
n_det = size(dets,1);
% [TODO] parsing only compatible with n_chrompohore = 2 for now
dets_char = parse_determinants(dets,num_chromophore); 

%visualize_CI_vec(dets,dets_char,3)

idx_S = find(dets_char == 1);
idx_CT = find(dets_char == 2);
idx_TT = find(dets_char == 3);

Proj = zeros(3,n_det);
Proj(1,idx_S) = 1/sqrt(length(idx_S));
Proj(2,idx_CT) = 1/sqrt(length(idx_CT));
Proj(3,idx_TT) = 1/sqrt(2);%1/sqrt(length(idx_TT)); 

J0 = -1000;
J_interact = [0, J0, 1/sqrt(12)*J0;
              J0, 0, J0/2;
              1/sqrt(12)*J0, -J0/2, 0];
J_coul = 5000;

%visualize_CI_vec(dets(idx_TT,:),dets_char(idx_TT),2)

%% Build FCI Hamiltonian
Hfci = build_fci_model_hamiltonian(dets,dets_char,J_interact,J_coul,omega_ge);
[CI_vec, D] = eig(Hfci);
E = diag(D);
subplot(2,2,[2,4])
plot(1:length(E),E,'bo')

%% Project onto S, CT, and TT subspace
PHP = Proj*Hfci*Proj';
[CI_vec_P, D] = eig(PHP);
EP = diag(D);
PHP
CI_vec_P;

for i = 1:3
    [~,idx(i)] = max(abs(CI_vec_P(i,:)));
end
% assign S (=1), CT (=2), and TT (=3) character to each eigenvector
CI_vec_P = CI_vec_P(:,idx);

EP = EP(idx)

subplot(221)
imagesc([1,2,3],[1,2,3],abs(CI_vec_P))
title(sprintf('E(S) = %4.2f, E(CT) = %4.2f, E(TT) = %4.2f',EP(1),EP(2),EP(3)));
colorbar
colormap(bluewhitered)
set(gca,'FontSize',14,'Linewidth',2,'xtick',[1,2,3],'xticklabels',{'S','CT','TT'},'ytick',[1,2,3],'yticklabels',{'S','CT','TT'})
CI_vec_P
    
%% Solve redfield equations
c = 3e-2; % cm/ps
kbT = 0;
npts_per = 20;
tau_max = 200; 
Nteval = tau_max*npts_per; dt = tau_max/Nteval;
tau = 0:dt:tau_max; % Ha^-1
taus = tau/c; % ps

omega_max = max(EP);
domega = omega_max/512;
omega_c = omega_max/25*ones(1,size(PHP,1));
gamma_i = 1*ones(1,size(PHP,1));

xin = [omega_c; gamma_i]; xin = xin(:);
fprintf('   Constructing Redfield tensor...\n')
tic
[Red_exc, GammaPop, JJ, Q] = RedfieldTensor(xin, domega, omega_max, tau_max, npts_per, EP, CI_vec_P, kbT);
fprintf('   Redfield tensor constructed in %4.1f seconds\n',toc)
subplot(222)
plot(domega:domega:omega_max,JJ{1},'Linewidth',2)
xlabel('\omega_{sys} / cm^{-1}')
ylabel('System-Bath Interaction Energy / cm^{-1}')
title('Spectral Density')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')

GammaPop = real(GammaPop)

%% Integrate equations and plot

rho = zeros(3,length(tau));
rho(:,1) = [1,0,0];
for i = 2:length(tau)
    rho(:,i) = expm(GammaPop*tau(i))*rho(:,1);
end

su
subplot(2,2,[3:4])
plot(taus,rho,'Linewidth',2)
xlabel('t / ps')
ylabel('Population')
title('Redfield Dynamics')
ll = legend('Singlet','CT','Triplet');
set(ll,'FontSize',13,'Location','NorthEast');
set(gca,'FontSize',14,'Linewidth',2,'Box','off')

figure(3)
for p = 1:3
    for q = 1:3
        for r = 1:3
            for s = 1:3
                plot(tau,real(squeeze(Q(p,q,r,s,:))))
                hold on
            end
        end
    end
end
hold off


