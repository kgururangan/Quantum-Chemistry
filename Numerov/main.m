clear all
clc
close all

%load 631gs_h2o.mat
load rhf_H2.mat
%
Ebind = Etot;
%RJ = theta_rad;

Hatowv =  219474.631370;  
bohrtoang = 0.529177;
[Ue,ie] = min(Ebind);
x = (RJ-RJ(ie))*bohrtoang; U = (Ebind + abs(Ue))*Hatowv; dx = abs(x(2)-x(1)); Nx = length(x);

mass_H = 1.00784; % amu
mA = mass_H; mB = mass_H; 
mu = (mA*mB)/(mA+mB);
%Elim = [10000, 55000];
Elim = [];
num_nodes = [0:2,4:20];
%num_nodes = 2;
flag_verbose = 1;
flag_plot = 0;
maxit = 200;
maxreit = 12;
tol = 1e-9;

[energy, psi, energy_hist, delta_hist, flag_exit] = numerov_wrap(x,U,mu,Elim,num_nodes,maxit,tol,maxreit,flag_verbose);
%[energy2, psi2, energy_hist, delta_hist, flag_exit] = numerov_function_v3(x,U,mu,Elim,num_nodes,maxit,tol,maxreit,flag_verbose);

%
%[energy, psi, energy_hist, delta_hist, flag_exit] = numerov_function_v1(x,U,mu,Elim,num_nodes,maxit,flag_verbose,flag_plot);

%%
plot(x,psi,x,psi2)
%%
iconv = find(flag_exit == 0);
plotWF(x,psi,U,energy,iconv);

%% Fit to spectroscopic parameters
morsefcn = @(x) x(1)*(iconv-1+1/2) + x(2)*(iconv-1+1/2).^2;

x0 = [5000,200]';
lb = []; ub = [];
fun = @(x) 1/length(iconv)*(energy(iconv)-morsefcn(x)).^2;
options = optimoptions('lsqnonlin','Display','Iter-Detailed');
[xsol,fval] = lsqnonlin(fun,x0,lb,ub,options);
%leg = cell(1,1)

% birge-sponer plot
Gv = []; ct = 1; nubs = [];
for i = 2:length(iconv)
    if iconv(i)-iconv(i-1) == 1
        Gv(ct) = energy(iconv(i))-energy(iconv(i-1));
        nubs(ct) = iconv(i-1);
        ct = ct + 1;
    end
end
c1 = polyfit(nubs,Gv,1);

figure(2)
subplot(121)
h1 = plot(iconv-1,morsefcn(xsol),'k-.','color',[0.8,0,0],'Linewidth',2); hold on
scatter(iconv-1,energy(iconv),100,'MarkerFaceColor',[0.8,0,0]); hold off
leg{1} = sprintf('\\omega_e = %4.2f cm^{-1}\n\\omega_ex_e = %4.2f cm^{-1}', xsol(1), -xsol(2));
ll = legend(h1,leg); set(ll,'FontSize',20,'Location','SouthEast');
xlabel('\nu (state)')
ylabel('E(\nu) / cm^{-1}')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')
axis tight
grid on
title('Energy vs. Vibrational State')

subplot(122)
scatter(nubs,Gv,100,'MarkerFaceColor',[0,0,0.8]); hold on
h2 = plot(nubs,c1(1)*nubs+c1(2),'k-.','color',[0,0,0.8],'Linewidth',2); hold off
leg{1} = sprintf('\\omega_ex_e = %4.2f cm^{-1}',-c1(1)/2);
ll = legend(h2,leg); set(ll,'FontSize',20,'Location','NorthEast');
xlabel('\nu (state)')
ylabel('G(\nu) = E(\nu+1)-E(\nu) / cm^{-1}')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')
axis tight
grid on
title('Birge-Sponer Plot')
