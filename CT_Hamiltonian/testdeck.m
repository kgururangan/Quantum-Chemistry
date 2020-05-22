clear all
clc
close all

% parameters
omega0 = 1000;
c = 3e-5;

par_CT.te = omega0;
par_CT.th = omega0;
par_CT.num_sites = 2;
par_CT.omega_ge = 12500; % singlet energy gap
par_CT.omega_ct = 12700;
par_CT.ct_range = par_CT.num_sites;

par_vib.n_quanta = 2;
par_vib.omega_g = [omega0];
par_vib.Dpv0 = [0.2];
par_vib.Dpv_pos = par_vib.Dpv0/2;
par_vib.Dpv_neg = par_vib.Dpv0/2;
par_vib.anharmonic_tensor = zeros(length(par_vib.omega_g),length(par_vib.omega_g),length(par_vib.omega_g));

par_spec.Nw = 1024;
par_spec.bw_Elec = 8000;
par_spec.ptol = 1e-6;
par_spec.kbT = 20;
par_spec.gamma = 0.28*omega0;
par_spec.cent_Elec = par_CT.omega_ge + omega0;
par_spec.popmin = 0.0;

% PROBLEM: In Frenkel and CT Hamiltonians, Hvib (which includes
% displacement) is only added for for an occupied site. By definition, the
% ground state is the fully unocciped vacuum which vibrations, so the 0
% quanta ground state cannot connect with any


par_CT.Jct = par_CT.te*par_CT.th/(par_CT.omega_ct - par_CT.omega_ge);
par_CT.Jct = 1.0*omega0;
par_CT.Jcoul = 0.0*omega0;

%states = get_states(par_CT.num_sites);

[I1, omega, MatEG, omega_mat, Ve, Vg, Ee, Eg, H_e, H_g] = run_CT_absorption(par_CT, par_vib, par_spec);

% par_dyn.omega_c = omega0*ones(1,size(H_e,1));
% par_dyn.gamma_i = 0.05*ones(1,size(H_e,1));
[taus, tau] = time_from_freq(omega, c);

% omega_c = omega0*ones(1,size(H_e,1));
% gamma_i = 0.05*ones(1,size(H_e,1));
% xin = [omega_c; gamma_i]; xin = xin(:);
% fprintf('   Constructing Redfield tensor...\n')
% tic
% [Red_exc, JJ, omega] = RedfieldTensor(xin, tau, Ee, Ve, par_spec.kbT);
% fprintf('   Redfield tensor constructed in %4.1f seconds\n',toc)

%%


%%
beta = linspace(0,3,5); 
Itot = cell(1,length(beta));
for i = 1:length(beta)
    par_CT.Jcoul = beta(i)*omega0;
    [I1, omega, MatEG, omega_mat, Ve, Vg, Ee, Eg, H_e, H_g] = run_CT_absorption(par_CT, par_vib, par_spec);
    Itot{i} = I1;
end

%%

close all

hold on
omegap = (omega - par_CT.omega_ge)/par_vib.omega_g(1);
for i = 1:length(Itot)
    plot(omegap,Itot{i})
end
set(gca,'FontSize',14,'Linewidth',2)
xlabel('(\omega - \omega_{ge})/\omega_{0g}')
hold off

%%