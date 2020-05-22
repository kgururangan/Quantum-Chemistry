clear all
clc
close all

%
addpath('/Users/karthik/Documents/MATLAB/2DPES/QPT');
addpath('/Users/karthik/Documents/MATLAB/Colormaps');
addpath('/Users/karthik/Documents/MATLAB/Additional_Functions');
addpath('/Users/karthik/Documents/MATLAB/Additional_Functions/regu');
%

cmap = jet_white;
pltc = default_cmap();

% problem dimension parameters
n_quanta = 2;
n_modes = 1;                    % need to fix energies for modes >1... are you adding omega_ge to the EC matrix too many times?
N = n_quanta^n_modes;           % dimension of vibrational basis per manifold
numstates = 3;                  % full vibronic basis dimension (product states of |elec>|vib>)
Nproblem = numstates*N;
% Nproblem = N^numstates;
zp = 0;
kbT = 20;
hbar = 1;
c = 3e-5;              % speed of light in cm/fs

%bandgap = 5000;
eegap = [800];
bandgap = 12500;
omega_g = [215];
omega_ge = [0, bandgap, bandgap + eegap];
omega_ge = omega_ge(1:numstates);
Bv = (1 + 1e-6)*ones(numstates,n_modes);
Bv(1,:) = ones(1,n_modes);


% axes
N2 = 64;
N3 = 128;
N4 = 64;

%bw_Elec = 2000;
bw_Elec = 1500;
bw_Raman = 1500;
dom2 = bw_Elec/N2;
dom3 = bw_Raman/N3;
dom4 = bw_Elec/N4;
cent_Elec = bandgap + eegap(1)/2;

omega_2 = (-dom2*(N2-1):2*dom2:dom2*(N2-1)) + cent_Elec;
omega_3 = -dom3*(N3-1):2*dom3:dom3*(N3-1);
omega_4 = (-dom4*(N4-1):2*dom4:dom4*(N4-1)) - cent_Elec;
tau_3 = 2*pi*(0:1/2/dom3/(N3):(N3-1)*1/2/dom3/N3); % [cm]
tau_4 = 2*pi*(0:1/2/dom4/(N4):(N4-1)*1/2/dom4/N4); % [cm]
tau_2 = 2*pi*(0:1/2/dom2/(N2):(N2-1)*1/2/dom2/N2); % [cm]
taus = tau_3/c;  % [fs]


tic
%%%%%%%%%%%%%%%%%%%%% Start of spectrum simulation %%%%%%%%%%%%%%%%%%%%%%%

% Step 1: define the SHO vibrational basis operators
Cre = zeros(n_quanta,n_quanta);
Ann = zeros(n_quanta,n_quanta);

for i = 0:n_quanta-1
    for j = 0:n_quanta-1
        
        if j == i + 1
            Ann(i+1,j+1) = sqrt(j);
        end
        
        if i == j + 1
            Cre(i+1,j+1) = sqrt(i);
        end
    end
end
Num = Cre*Ann;

% Step 2: Ground state (0-exciton manifold)          
H0_g = zeros(N,N);
for i = 1:n_modes
    if zp ~= 0
        H0temp = omega_g(i)*(Num + 0.5*eye(n_quanta));
    else
        H0temp = omega_g(i)*Num;
    end
    A = 1;
    for j = 1:n_modes  
        if j ~= i
            A = kron(eye(n_quanta),A);
        else
            A = kron(H0temp,A);
        end
    end
    H0_g = H0_g + A;
end
[V_g,E_g] = eig(H0_g);

% Generate Redfield data
M1 = 40;
M2 = 40;
MJ = 5;
MG = 60;

%displ = linspace(0,1.25,M); displ1 = displ; displ2 = displ;
displ1 = [0, 0.1, 0.5, 1.2];
displ2 = [0, 0.1, 0,5, 1.2];
coulcoupling = linspace(20,200,MJ);
vibenergy = [linspace(80,400,20),linspace(401,1800,MG),linspace(1801,2000,20)];

[~,A1fix] = min( abs(displ1 - 1.2) );
[~,A2fix] = min( abs(displ2 - 0.1) );
[~,A3fix] = min( abs(coulcoupling - 100) );
[~,A4fix] = min( abs(vibenergy - 215) );

A1par = 1:length(displ1);
A2par = 1:length(displ2);
A3par = A3fix;
A4par = 1:length(vibenergy);

M1 = length(A1par); M2 = length(A2par); MJ = length(A3par); MG = length(A4par);

RED_EXC = cell(M1,M2,MJ,MG);
GammaPop = cell(M1,M2,MJ,MG);
H_E = cell(M1,M2,MJ,MG);
XP = cell(M1,M2,MJ,MG);

for A1 = A1par
    for A2 = A2par
        for A3 = A3par
            for A4 = A4par
        
                Dpv = [0;displ1(A1);displ2(A2)];
                
                Jcoul = zeros(numstates-1);
                Jcoul(1,2) = coulcoupling(A3); Jcoul(2,1) = coulcoupling(A3);
                omega_e = vibenergy(A4);
                
                % Step 3: Excited states
                H_e = zeros(N*(numstates-1));
                Ue = cell(1,numstates-1);
                for k = 2:numstates
                    for kp = 2:numstates
                        if k == kp
                            Htemp = zeros(N);
                            for i = 1:n_modes
                                %Ann_e = (Ann - Dpv(k,i)*eye(n_quanta));
                                %Cre_e = Ann_e';
                                %Num_e = Cre_e*Ann_e;
                                if zp ~= 0
                                    H0temp = omega_e(i)*(Num_e + 0.5*eye(n_quanta));
                                else
                                    %H0temp = Bv(k,i)^2*omega_g(i)*Num_e;
                                    H0temp = omega_e(i)*Num+omega_e(i)*Dpv(k,i)*(Ann+Cre);
                                end
                                A = 1;
                                for j = 1:n_modes
                                    if j ~= i
                                        A = kron(eye(n_quanta),A);
                                    else
                                        A = kron(H0temp,A);
                                    end
                                end
                                Htemp = Htemp + A;
                            end
                            Htemp = Htemp + omega_ge(k)*eye(N);
                            kernel  = zeros(numstates-1); kernel(k-1,k-1) = 1;
                            [vtemp,etemp] = eig(Htemp);
                            Ue{k-1} = vtemp;
                            %H_e = H_e + kron(kernel,etemp);
                            H_e = H_e + kron(kernel,Htemp);
                        else
                            Htemp = Jcoul(k-1,kp-1)*eye(N);
                            kernel = zeros(numstates-1); kernel(k-1,kp-1) = 1;
                            H_e = H_e + kron(kernel,Htemp);
                        end
                    end
                end
                
                Ne = length(H_e);
                
                [V_e,D_e] = eig(H_e);
                
                omega_c = 200*ones(1,Ne);
                gamma_i = 0.05*ones(1,Ne);
                xin = [omega_c; gamma_i]; xin = xin(:);
                tic
                [Red_exc, JJ, omega] = RedfieldTensor(xin, tau_3, diag(D_e), V_e, kbT);
                toc
                
                % Store the results
                XP{A1,A2,A3,A4} = [displ1(A1), displ2(A2), coulcoupling(A3), vibenergy(A4), Ne, omega_c, gamma_i];
                RED_EXC{A1,A2,A3,A4} = Red_exc;
                H_E{A1,A2,A3,A4} = H_e;
                
                R = reshape(Red_exc,Ne^2,Ne^2);
                S = zeros(Ne);
                for p = 1:Ne
                    for q = 1:Ne
                        I1 = ind2D(Ne,p,p);
                        I2 = ind2D(Ne,q,q);
                        S(p,q) = R(I1,I2);
                    end
                end
                GammaPop{A1,A2,A3,A4} = real(S);
            end
        end
    end
end

%% Rate contours
Ne = length(H_e);
rates = zeros(Ne,Ne,M1,M2,MJ,MG);
for i = A1par
    for j = A2par
        for m = A3par
            for n = A4par
                
                R = GammaPop{i,j,m,n};
                for k = 1:Ne
                    for l = 1:Ne
                        rates(k,l,i,j,m,n) = R(k,l);
                    end
                end
                
            end
        end
    end
end

%parameters for figure and panel size
plotheight=100;
plotwidth=100;
subplotsx=4;
subplotsy=4;   
leftedge=2;
rightedge=2;   
topedge=2;
bottomedge=2;
spacex=0.3;
spacey=5;
fontsize=5;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
 
%setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
 
%loop to create axes
for i=1:subplotsx
    for ii=1:subplotsy
        
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        contourf(displ2,displ1,real(squeeze(rates(i,ii,:,:,A3fix,A4fix))));  hold on;
        cc =colorbar; set(cc,'Location','eastoutside');
        cc.FontSize = 12;
        ll = legend(sprintf('%d -> %d',ii,i)); set(ll,'FontSize',13,'Location','NorthEast');
        set(gca,'FontSize',13,'Ydir','normal','Linewidth',2,'Box','off')
        colormap(jet_white)
        axis square
        grid on
        
        xlabel('d_2'); ylabel('d_1');
        
    end
end

%% Find rhoM using Redfield tensor

rho33pop0 = 0.50;
MI = length(rho33pop0);

RHO = cell(M1,M2,MJ,MG,MI);
TAU_RISE = cell(1,MI);
C0 = cell(1,MI);
C1 = cell(1,MI);
FIT = cell(1,MI);
GOPT = cell(1,MI);

for init = 1:MI

rho_site0 = zeros(Ne); 
rho_site0(Ne-1,Ne-1) = rho33pop0(init); rho_site0(Ne,Ne) = 1-rho_site0(Ne-1,Ne-1);

tic
for i = A1par
    for j = A2par
        for l = A3par
            for m = A4par
      
                H_e = H_E{i,j,l,m};
                [V_e,D_e] = eig(H_e);
                % Integrate in time
                Nteval = 100;
                Ne = length(H_e);
                teval = linspace(min(tau_3),max(tau_3),Nteval);
                R_rs = reshape(Red_exc,Ne^2,Ne^2);
                rho0 = V_e'*rho_site0*V_e;
                rhoM = zeros(Ne^2,Nteval);
                for k = 1:length(teval)
                    rhoM(:,k) = expm(R_rs*teval(k))*rho0(:);
                end 
                % Pluck out population evolution
                rho_pop = zeros(Ne,Nteval);
                for k = 1:Ne
                    rho_pop(k,:) = real(rhoM(ind2D(Ne,k,k),:));
                end
                RHO{i,j,l,m,init} = rho_pop;
                
            end
        end
    end
end
toc
        
% Fit rho_11 evolution
options = optimoptions('LSQNONLIN','Display','iter-detailed','Diagnostics','on');
options.FunctionTolerance = 1e-15;
options.OptimalityTolerance = 1e-20;
options.StepTolerance = 1e-20;

tau_rise = zeros(M1,M2,MJ,MG);
tau_fall = zeros(M1,M2,MJ,MG);
c1 = zeros(M1,M2,MJ,MG);
c2 = zeros(M1,M2,MJ,MG);
c0 = zeros(M1,M2,MJ,MG);
fit = cell(M1,M2,MJ,MG);
gopt = zeros(M1,M2,MJ,MG);

biexp = 0; 

for i = A1par
    for j = A2par
        for l = A3par
            for m = A4par
        
                rho11 = RHO{i,j,l,m}(1,:);
                
                if biexp == 1 % c0 + c1*exp(-a*t) - c2*exp(-b*t)
                    nparam = 5;
                    model = @(x) x(1) + x(2)*exp(-x(3)*teval) - x(4)*exp(-x(5)*teval);
                    fun = @(x) 1/Nteval*abs(model(x) - rho11) ;
                    lb = zeros(1,nparam);
                    ub = [1,10,50,10,50];
                    x0 = (ub-lb).*rand(1,nparam) + lb;
                    [xsolve,fval] = lsqnonlin(fun,x0,lb,ub,options);
                    gopt(i,j,l,m) = fval;
                    fit{i,j,l,m} = model(xsolve);
                    tau_rise(i,j,l,m) = xsolve(3); tau_fall(i,j,l,m) = xsolve(5);
                    c1(i,j,l,m) = xsolve(2); c2(i,j,l,m) = xsolve(4); c0(i,j,l,m) = xsolve(1);
                else % c0 + c1*(1 - exp(-a*t))
                    nparam = 3;
                    model = @(x) x(1)*(1-exp(-x(2)*teval))+x(3);
                    %model = @(x) rho11(1) + rho11(end)*(1-exp(-x*teval));
                    fun = @(x) 1/Nteval*abs(model(x) - rho11) ;
                    lb = zeros(1,nparam);
                    ub = [1,50,1];
                    %ub = 50;
                    x0 = (ub-lb).*rand(1,nparam) + lb;
                    [xsolve,fval] = lsqnonlin(fun,x0,lb,ub,options);
                    gopt(i,j,l,m) = fval;
                    fit{i,j,l,m} = model(xsolve);
                    %tau_rise(i,j,l,m) = xsolve;
                    tau_rise(i,j,l,m) = xsolve(2); c1(i,j,l,m) = xsolve(1); c0(i,j,l,m) = xsolve(3);
                end
            
            end
        end
    end
end

TAU_RISE{init} = tau_rise; C0{init} = c0; C1{init} = c1; FIT{init} = fit; GOPT{init} = gopt;

end
%% Check terms in fit
clf
close all
k1 = 1; k2 = 1; k3 = A3fix; k4 = 10;
plot(teval,RHO{k1,k2,k3,k4,1}(1,:),teval,fit{k1,k2,k3,k4},'Linewidth',2)
grid on
ll = legend('Data','Fit');
%% Displacement Contour Plot
close all
%tau_rise_plot = TAU_RISE{2};
%tau_rise_plot = imgaussfilt(tau_rise,1);
[a,b] = contourf(displ,displ,squeeze(tau_rise(:,:,A3fix,A4fix)),1e3); 
set(b,'EdgeColor','None')
%imagesc(displ,displ,tau_rise_plot); hold on;
%contour(displ,displ,tau_rise_plot,15,'b-','Linewidth',2)
axis square
grid on
set(gca,'FontSize',20,'Ydir','normal','Linewidth',2,'Box','off')
colorbar
colormap(jet_white)
xlabel('d_2')
ylabel('d_1')
title('Initial Condition: \rho_{33}(0) = 0.5, \rho_{44}(0) = 0.5')

%% J-coupling curve
plot(coulcoupling,squeeze(tau_rise(A1fix,A2fix,:,A4fix)))
%% omega_e curve
%clr = brewermap(8,'Set3');
clr = jet(8);
leg = cell(6,1);
hold on;
count = 1;
for i = [1,2,4]
    for j = [1,2]
        for k = A3par
                plot(vibenergy/(eegap+coulcoupling(A3fix)),squeeze(tau_rise(i,j,k,:)),'bo','MarkerSize',8,'MarkerFaceColor',clr(count,:))
                %plot(vibenergy/(eegap+2*coulcoupling(A3fix)),squeeze(tau_rise(i,j,k,:)),'k-.','color',clr(count,:))
                leg{count} = ['d_1 = ',num2str(displ1(i)),' , ','d_2 = ',num2str(displ2(j))];
                count = count + 1;
        end
    end
end
grid on
set(gca,'FontSize',17,'Linewidth',2,'Box','off')
xlabel('\omega_{e}/(\Delta_{ee}+2J)')
ylabel('k_{tr} [cm^{-1}]')
title('Initial Condition: \rho_{33}(0) = 0.5, \rho_{44}(0) = 0.5')
%ll = legend(['J = 10 cm^{-1}'],['J = 57.5 cm^{-1}'],['J = 105 cm^{-1}'],['J = 152.5 cm^{-1}'],['J = 200 cm^{-1}']); 
%ll = legend(['d_1 = 1.2' newline 'd_2 = 0.05' newline 'J = 10 cm^{-1}'],['d_1 = 1.2' newline 'd_2 = 0.05' newline 'J = 200 cm^{-1}']); 
ll = legend(leg);
set(ll,'Fontsize',16,'Autoupdate','off');
%%
count = 1;
for i = A1par
    for j = A2par
        for k = 1:MJ
                %plot(vibenergy/eegap,squeeze(tau_rise(i,j,k,:)),'bo','MarkerSize',8,'MarkerFaceColor',clr(count,:))
                plot(vibenergy/eegap,squeeze(tau_rise(i,j,k,:)),'k-.','color',clr(count,:),'Linewidth',2)
                count = count + 1;
        end
    end
end
%% d1/d2 vs. omega_e/Delta contour
close all
tau_rise_plot = imgaussfilt(tau_rise,1);
[a,b] = contourf(vibenergy(:)/eegap,displ,squeeze(tau_rise(:,A2fix,A3fix,:)),1e3); 
set(b,'EdgeColor','None')
axis square
grid on
set(gca,'FontSize',20,'Ydir','normal','Linewidth',2,'Box','off')
colorbar
colormap(jet)
xlabel('\omega_e/\Delta_{ee}')
ylabel('d_1/d_2')
title('Initial Condition: \rho_{33}(0) = 0.8, \rho_{44}(0) = 0.2')