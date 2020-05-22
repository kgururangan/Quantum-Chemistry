clear all
clc
close all

modes = [215];
numspec = length(modes);
Si = cell(numspec,1);

for NUM = 1:numspec

tic

flag.ESE = 1; 
flag.GSB = 0; 
flag.ESA = 0;  
redfield = 0; % only applies to ESE for now

%
% addpath('/Users/karthik/Documents/MATLAB/2DPES/QPT');
% addpath('/Users/karthik/Documents/MATLAB/Colormaps');
% addpath('/Users/karthik/Documents/MATLAB/Additional_Functions');
% addpath('/Users/karthik/Documents/MATLAB/Additional_Functions/regu');
%
cmap = parula;
pltc = jet;

% problem dimension parameters
n_quanta = 3;
n_modes = 1;                    % need to fix energies for modes >1... are you adding omega_ge to the EC matrix too many times?
N = n_quanta^n_modes;           % dimension of vibrational basis per manifold
numstates = 3;                  % full vibronic basis dimension (product states of |elec>|vib>)
Nproblem = numstates*N;
zp = 0;
hbar = 1;
c = 3e-5;              % speed of light in cm/fs
                       
% problem parameters
Jcoul = zeros(numstates); % vibronic coupling constant, sign doesn't seem to matter [eV]
indtr_state = [4,7]; % decide which states transfer via J-coupling HERE
% Jcoul < 0 = downhill transfer, > 0 = uphill transfer
Jcoul(2,3) = 200; Jcoul(3,2) = 200;
% Jcoul = zeros(numstates)

%bandgap = 5000;
eegap = [800];
bandgap = 12500;
%omega_g = [215];
omega_g = modes(NUM);
omega_ge = [0, bandgap, bandgap + eegap];
Bv = (1 + 1e-6)*ones(numstates,n_modes);
Dpv = zeros(numstates,n_modes);
% first row of Bv and Dpv MUST be 0's to get the correct ground state block!
Bv(1,:) = ones(1,n_modes);
Dpv(1,:) = zeros(1,n_modes);
Dpv = [0; 0.2; 0.2]; % changes uphill and downhill??

% spectrum parameters
pmin = 1e-4;
pop_min = 0.5;
kbT = 20;

% axes
N2 = 100;
N3 = 300;
N4 = 100;

%bw_Elec = 2000;
bw_Elec = 1000;
bw_Raman = 1000;
dom2 = bw_Elec/N2;
dom3 = bw_Raman/N3;
dom4 = bw_Elec/N4;
cent_Elec = bandgap + eegap(1)/2;

omega_2 = (-dom2*(N2-1):2*dom2:dom2*(N2-1)) + cent_Elec;
omega_3 = -dom3*(N3-1):2*dom3:dom3*(N3-1);
omega_4 = (-dom4*(N4-1):2*dom4:dom4*(N4-1)) - cent_Elec;
tau_3 = 2*pi*(0:1/2/dom3/(N3):(N3-1)*1/2/dom3/N3); % [cm]
taus = tau_3/c;  % [fs]


% lineshape parameters - discrete lifetimes
Gamma_sqc_ge2 = 100;     % Single quantum coherence (between g and e2)
Gamma_sqc_ge1 = 100;     % Single quantum coherence (between g and e1)
%Gamma_zqc_e1e2 = 70;    % Electronic coherence ZQC e1 and e2, sets lifetime of electronic crosspeaks
%Gamma_zqc = 10;          % Zero quantum coherence

% lineshape parameters - continuous lifetimes
lifetime_SQC_0 = 25; % lifetime of SQC state between g and e1 [fs]
lifetime_POP_0 = 1e6;  % lifetime of population state [fs]
Gamma_SQC_0 = (2*pi*c*lifetime_SQC_0)^-1; % lifetime of single quantum coherence between g and e1
Gamma_POP_0 = (2*pi*lifetime_POP_0)^-1;  % 
GAMMA = Gamma_SQC_0*bandgap;

par.tau_3 = tau_3;
par.omega_g = omega_g;
par.omega_2 = omega_2; par.N2 = N2;
par.omega_3 = omega_3; par.N3 = N3;
par.omega_4 = omega_4; par.N4 = N4;
par.bw_Elec = bw_Elec;
par.cent_Elec = cent_Elec;
par.bw_Raman = bw_Raman;
par.n_quanta = n_quanta; par.n_modes = n_modes; par.numstates = numstates;
par.Dpv = Dpv; par.Bv = Bv;
par.Gamma_sqc_ge1 = Gamma_sqc_ge1; par.Gamma_sqc_ge2 = Gamma_sqc_ge2;
par.omega_g = omega_g; par.omega_ge = omega_ge; par.eegap = eegap; par.bandgap = bandgap;
par.Jcoul = Jcoul;
par.kbT = kbT; par.pop_min = pop_min; par.pmin = pmin;
par.hbar = hbar; par.c = c;

% 'lorentzian', 'gaussian', 'voigt'
lineshape = 'gaussian';


tic
%%%%%%%%%%%%%%%%%%%%% Start of spectrum simulation %%%%%%%%%%%%%%%%%%%%%%%



Htot = zeros(Nproblem,Nproblem);
%FC_ge = cell(nelecstates,nelecstates);
FCmat = zeros(Nproblem,Nproblem);

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

% Step 2: Full vibronic Hamiltonian matrix

% Create vibronic Hamiltonian matrix by looping thorough all electronic
% states, constructing each as a NxN block in the vibrational basis
for k = 1:numstates
    for kp = 1:numstates
          
        if k == kp  % ground state SHO and excited state DHO blocks
            
            H = zeros(N,N);
            
            for i = 1:n_modes
                if zp ~= 0
                    H0temp = omega_g(i)*(Num + 0.5*eye(n_quanta));
                else
                    % effect of displacement on excited states taken care
                    % of by FC matrix
                    %H0temp =  omega_ge(k)*eye(n_quanta) +...
                    %    Bv(k,i)^2*omega_g(i)*( (Num + Dpv(k,i)^2/2*eye(n_quanta)) - (Cre + Ann)/sqrt(2)*Dpv(k,i) );
                    H0temp = Bv(k,i)^2*omega_g(i)*Num;
                end
                A = 1;
                for j = 1:n_modes
                    
                    if j ~= i
                        A = kron(eye(n_quanta),A);
                    else
                        A = kron(H0temp,A);
                    end
                    
                end
                H = H + A;
            end
            
            % no FC transition matrix for GG,E1E1,E2E2, etc.
            GmnE = zeros(N,N);
            %GmnE = myFranckCondon((1+1e-6)*ones(1,n_modes),Dpv(k,:),n_modes,n_quanta);
            H = H + omega_ge(k)*eye(N);

            %FC_ge{k,kp} = zeros(N,N);
            
            %if any(isnan(H))
            %    crash = [k,kp];
            %    break;
            %end
            
            
        else  
          
        % off-diagonal blocks, proportional to FC matrix
        % effecively an overlap approximation to non-adiabatic tensor
            rel_disp = abs(Dpv(k,:) - Dpv(kp,:));
            rel_beta = Bv(k,:)./Bv(kp,:) + 1e-6;
            GmnE = myFranckCondon(rel_beta, rel_disp, n_modes, n_quanta);
            %H = Jcoul(k,kp)*GmnE;
            H = Jcoul(k,kp)*eye(N);
            %FC_ge{k,kp} = GmnE;
            
            %if any(isnan(H))
            %    crash = [k,kp];
            %    break;
            %end
            
        end
        
           
        % kron into full vibronic H matrix and full FC matrix
        kernel = zeros(numstates);
        kernel(k,kp) = 1;
        Htot = Htot + kron(kernel,H);
        FCmat = FCmat + kron(kernel,GmnE);
     
    end
        
end

% obtain the diabatic basis functions and energies
[V,D] = eig(Htot);
Enn = diag(D);
[Ennsort,Isort] = sort(Enn);
Vsort = V(:,Isort);
Dsort = D(Isort,Isort);
breakdown = Vsort.^2;

% Step 3: Thermal Density Matrix
Z = trace(expm(-Dsort/kbT)); 
rho_eq = expm(-Dsort/kbT)/Z;
drho_eq = zeros(Nproblem,Nproblem);
for i = 1:Nproblem
    for j = 1:Nproblem   
        drho_eq(i,j) = (1/Z)*(exp(-Dsort(i,i)/kbT) - exp(-Dsort(j,j)/kbT));
    end
end

% Step 4: Transition Matrix Elements for ORF
% rotate from adiabatic to diabatic basis
%!!! SO MANY OF YOUR PROBLEMS CAN BE SOLVED HERE !!!
% <A|B> = sum_{G,E}(<A|G><G|E><E|B>) = Vsort'*FCmat*Vsort
MatGE = FCmat*Vsort; 
MatEG = MatGE';


% Step 5: Optical Response Functions
% alpha states
[Ialph,Jalph] = find(rho_eq >= pop_min); 
alpha = 1;
% ground states
[Ig] = find(Ennsort <= bw_Raman);     
% 1-exc states
[Ie] = find(Ennsort <= cent_Elec + bw_Elec/2 & Ennsort >= cent_Elec - bw_Elec/2); 
% f-states.. what's a good way to get f-states?
[~,Ifmin] = min( abs( Ennsort - abs((4*bw_Elec - cent_Elec)) ) );
[~,Ifmax ] = min( abs( Ennsort - abs((4*bw_Elec + cent_Elec))) );
Ifmin = Ie(end)+1;
If  = Ifmin:Ifmax;

E_ec = Ennsort(Ie); E_g = Ennsort(Ig); Ef = Ennsort(If);
Ne = length(Ie);    Ng = length(Ig);   Nf = length(If);  Nalph = length(Ialph); 

% Compute statics GM matrices
if flag.ESE == 1

    GMprocc = zeros(N2*N4,Ng,Ne,Ne,Ne,Ne);
    Ptese = zeros(Ng,Ne,Ne,Ne,Ne);
    
    for beta = 1:Ng
        for a = 1:Ne
            for b = 1:Ne
                for c = 1:Ne
                    for d = 1:Ne
                        
                        %pre = MatGE(1,Ie(a))*MatEG(Ie(c),Ig(beta))...
                        %*MatGE(Ig(beta),Ie(d))*MatEG(Ie(b),1)*rho_eq(1,1);
                        
                        pre = MatEG(Ie(a),1)*MatEG(Ie(b),1)*MatGE(Ig(beta),Ie(c))*MatGE(Ig(beta),Ie(d))*rho_eq(1,1);
                        %pre = 1;
                        Ptese(beta,a,b,c,d) = pre;
                       
                        % omega_2 spectrum
                        f2 = myLineshape_v3(omega_2,Ennsort(Ie(a))-Ennsort(1),Gamma_sqc_ge1,lineshape);
                        % omega_4 spectrum
                        f4 = myLineshape_v3(omega_4,Ennsort(Ig(beta))-Ennsort(Ie(d)),Gamma_sqc_ge1,lineshape);
                        
                        GMprocc(:,beta,a,b,c,d) = pre*kron(f4,f2);

                    end
                end
            end
        end
    end
    
    GMproc = reshape(squeeze(sum(GMprocc,2)),N2*N4,Ne^4);
    Ptese  = reshape(squeeze(sum(Ptese,1)),1,Ne^4);
    %Lkvib = reshape(squeeze(sum(sum(sum(GMprocc,3),5),6)),N2*N4,Ne^2);
    
    GMprocESE = GMproc;
    
else
    GMprocESE = zeros(N2*N4,Ne^4);
end

if flag.GSB == 1

    GMprocc = zeros(N2*N4,Ne,Ne,Ng,Ng,Ng,Ng);
    
    for a = 1:Ne
        for b = 1:Ne
            for eta = 1:Ng  % index running over alpha states
                for beta = 1:Ng
                    for gamma = 1:Ng
                        for delta = 1:Ng
                            
                            pre = MatGE(1,Ie(a))*MatGE(Ig(delta),Ie(b))...
                                *MatGE(Ig(beta),Ie(a))*MatEG(Ie(b),Ig(gamma))*rho_eq(1,1);
                            
                            % omega_2 spectrum
                            f2 = myLineshape_v3(omega_2,Ennsort(Ie(a))-Ennsort(1),Gamma_sqc_ge1,lineshape);
                            
                            %  if Ie(b) < Ie(ind_e2(1))
                            f4 = myLineshape_v3(omega_4,Ennsort(Ig(delta))-Ennsort(Ie(b)),Gamma_sqc_ge1,lineshape);
                            
                            GMprocc(:,a,b,eta,beta,gamma,delta) = pre*kron(f4,f2);
                            
                            
                        end
                    end
                end
            end
        end
    end
    
    GMproc = reshape(squeeze(sum(sum(GMprocc,2),3)),N2*N4,Ng^4);
    GMprocGSB = GMproc;
    
else
    GMprocGSB = zeros(N2*N4,Ng^4);
end

if flag.ESA == 1

    GMprocc = zeros(N2*N4,Ng,Nf,Ne,Ne,Ne,Ne);
    
    for beta = 1:Ng
        for f = 1:Nf
            for a = 1:Ne
                for b = 1:Ne
                    for c = 1:Ne
                        for d = 1:Ne
                            
                            pre = MatGE(1,Ie(a))*MatEG(Ie(b),1)...
                                *MatEG(If(f),Ie(d))*MatGE(Ie(c),If(f))*rho_eq(1,1);
                            
                            % omega_2 spectrum
                            f2 = myLineshape_v3(omega_2,Ennsort(Ie(a))-Ennsort(1),Gamma_sqc_ge1,lineshape);
                            
                            %  if Ie(b) < Ie(ind_e2(1))
                            f4 = myLineshape_v3(omega_4,Ennsort(Ie(c))-Ennsort(If(f)),Gamma_sqc_ge1,lineshape);
                            
                            GMprocc(:,beta,f,a,b,c,d) = pre*kron(f4,f2);
                            
                            
                        end
                    end
                end
            end
        end
    end
    
    GMproc = reshape(squeeze(sum(sum(GMprocc,2),3)),N2*N4,Ne^4);
    GMprocESA = GMproc;
    
else
    GMprocESA = zeros(N2*N4,Ne^4);
end

GM_toc = toc

if redfield == 1
    % Redfield dynamics (takes a long time..)
    
    omega_c = [10000];
    gamma_i = [0.02];
    %ind = [1,3];
    
    %omega_c = [100, 100, 400, 400];
    omega_c = repmat(omega_c,1,Ne);
    %gamma_i = [0.02, 0.02, 0.02, 0.02];
    gamma_i = repmat(gamma_i,1,Ne);
    
    xtrans = [omega_c;gamma_i];
    
    Ha = Htot(Ie,Ie);
    V_S = Vsort(Ie,Ie);
    E_S = Ennsort(Ie);
    VSS = kron(V_S,V_S);
    VSST = kron(V_S',V_S');
    
    tic
    [Red_exc, JJ, omega_bath] = RedfieldTensor(xtrans(:), tau_3, E_S, V_S, kbT);
    Redfield_toc = toc
    
    Red_exc_rs = reshape(Red_exc,Ne^2,Ne^2);
    Red_site_rs = VSS*Red_exc_rs*VSST;
    [TM,Ev] = TM_matrix_v2(Red_exc_rs,tau_3,0,'expm'); % use the excitonic redfield tensor for spectrum dynamics
    [TM_site,Ev_site] = TM_matrix_v2(Red_site_rs,tau_3,0,'expm');
    
    % GSB dynamics is matrix exponential still
    Gamma_zqc_e1e2 = 70;     % Electronic coherences
    Gamma_zqc = 10;          % Vibrational coherences
    Gammadab = [Gamma_zqc,Gamma_zqc_e1e2];
    Relax_gnd = diag(-Gammadab(1)*ones(1,Ng));
    GammaPop_gnd = zeros(Ng);
    Kpop_gnd = [0,  0;
                15,-15];
    indpos_gnd = [1,2];
    GammaPop_gnd(indpos_gnd,indpos_gnd) = Kpop_gnd;
    GammaPop_gnd = GammaPop_gnd + Relax_gnd;
    ind_0exc_man = ones(1,length(Ig));
    [TMGSB] = myDynamics_v5(par,Ennsort(Ig),[],GammaPop_gnd,Gammadab,ind_0exc_man);
else
    tic
    % Matrix exponential (phenomenological) dynamics
    
    Gamma_zqc_e1e2 = 10;     % Electronic coherences
    Gamma_zqc = 10;          % Vibrational coherences
    
    Gammadab = [Gamma_zqc,Gamma_zqc_e1e2];
    Relax_exc = diag(-Gammadab(1)*ones(1,Ne));
    Relax_gnd = diag(-Gammadab(1)*ones(1,Ng));
    
    GammaPop_exc = zeros(Ne);
    GammaPop_gnd = zeros(Ng);
    
    Kpop_exc = [0,  30;
                0,  -50];
    
   % Kpop_gnd = [0,  0;
     %           0,  0];
    
    indpos_exc = indtr_state - length(Ig);
   % indpos_gnd = [1,5];
    
    GammaPop_exc(indpos_exc,indpos_exc) = Kpop_exc;
    %GammaPop_gnd(indpos_gnd,indpos_gnd) = Kpop_gnd;
    
    GammaPop_exc = GammaPop_exc + Relax_exc;
    GammaPop_gnd = GammaPop_gnd + Relax_gnd;
    
    ind_1exc_man = ones(1,length(Ie));
    ind_0exc_man = ones(1,length(Ig));
    
    [TM] = myDynamics_v5(par,Ennsort(Ie),[],GammaPop_exc,Gammadab,ind_1exc_man);
    [TMGSB] = myDynamics_v5(par,Ennsort(Ig),[],GammaPop_gnd,Gammadab,ind_0exc_man);
    
    TM_toc = toc
end

% Add dynamics using GMproc
posv = zeros(Ne,Ne);
for a = 1:Ne
    for c = 1:Ne
        b = a; d = c;
        J = ind2D(Ne,a,b);
        K = ind2D(Ne,c,d);
        posv(c,a) = ind2D(Ne^2,J,K);
    end
end
pos = posv(:); % vector containing population terms in 1-exciton manifold

posvgnd = zeros(Ng,Ng);
for alpha = 1:Ng
    for gamma = 1:Ng
        beta = alpha; delta = gamma;
        J = ind2D(Ng,alpha,beta);
        K = ind2D(Ng,gamma,delta);
        posvgnd(gamma,alpha) = ind2D(Ng^2,J,K);
    end
end
posgnd = posvgnd(:); % vector containing population terms in GSM 

ESE = reshape(GMprocESE*TM,N2,N4,N3); ESEdyn = reshape(GMprocESE(:,pos)*TM(pos,:),N2,N4,N3);
GSB = reshape(GMprocGSB*TMGSB,N2,N4,N3); GSBdyn = reshape(GMprocGSB(:,posgnd)*TMGSB(posgnd,:),N2,N4,N3);
ESA = reshape(GMprocESA*TM,N2,N4,N3); ESAdyn = reshape(GMprocESA(:,pos)*TM(pos,:),N2,N4,N3);

I4WMt = ESE + GSB - ESA;
I4WMdyn = ESEdyn + GSBdyn - ESAdyn;

posDiff = zeros(1,Ne^2);
count = 1;
for i = 1:Ne
    for j = 1:Ne
        posDiff(count) = Ennsort(Ie(i))-Ennsort(Ie(j));
        count = count+1;
    end
end

Si{NUM} = I4WMt;

end



%%%%%%%%%%%%%%%%%%%%% End of spectrum simulation %%%%%%%%%%%%%%%%%%%%%%%

%%
coeffs = [1,0,0];
Stot = coeffs(1)*Si{1} + coeffs(2)*Si{2} + coeffs(3)*Si{3};

%% Look at the full 3D spectrum.. just for fun
vis_flag = 1;
I4WMt_fft = fft(I4WMt,N3,3);
cm = [0.3,0.4,0.8];
if vis_flag ~= 0
    [X,Y,Z] = meshgrid(omega_2,omega_4,taus);
    clf; hold on; view(3);
    C = redblue(1350558);
    %C = brewermap((205551),'RdBu');
    %C = flipud(jet(89865));
    %C = (cmocean('balance',1202202));
    %C = fire(1266006);
    %C = jet(1266006);
    p = patch(isosurface(X,Y,Z,abs(I4WMt)./max(max(max(abs(I4WMt)))),'noshare'),'FaceVertexCData',C);
    length(p.Vertices) % check this to set size of C
    p.FaceColor = 'none';
    p.EdgeColor = 'interp';
    view(50,40)
    %contourslice(X,Y,Z,abs(I4WMt),[],[],taus(1:20:N3),20); 
    colormap(jet_white)
    camlight
    %lighting gouraud
    grid('on')
    set(gca,'fontsize',24)
    axis('tight')
    xlabel('\omega_1 (cm^{-1})'); ylabel('\omega_3 (cm^{-1})'), zlabel('\tau_2 (fs)')
end

%%
close all

Stot = Si{1};

% answers blow up if tau_3 is too long? -> instability caused when unequal
% excited state displacements..

rho0_site = zeros(Ne,Ne);
% populationsx
rho0_site(1,1) = 0; rho0_site(2,2) = 0; rho0_site(3,3) = 0; rho0_site(4,4) = 1;
% coherences
rho0_site(1,2) = 0; rho0_site(2,1) = 0;

SPEC = I4WMdyn;

flag.plot_spectrum = 1;
flag.plot_trace = 0;
flag.scroll_spectrum = 0;
flag.rhomat = 0;
flag.plotsignals = 0;
flag.plot_time_trace = 0;


%v = VideoWriter('DAS.avi');
%v.FrameRate = 7;
%open(v);
if flag.plot_spectrum == 1
figure(3)
for i = 1:N3
    imagesc(-omega_4,omega_2,squeeze(real(SPEC(:,:,i)))); hold on;
    contourf(-omega_4,omega_2,squeeze(real(SPEC(:,:,i))),10);
    axis square
    title(sprintf('T = %4.2f ps', tau_3(i)/3e-2))
    xlabel('-\omega_t/cm^{-1}');
    ylabel('\omega_\tau/cm^{-1}');
    set(gca,'FontSize',20,'YDir','normal','Linewidth',1.5,'Box','off')
    %set(gcf,'color','w');
    colorbar
    %colormap(brewermap([],'RdBu'))
    colormap(cmap)
    caxis( [min(min(real(squeeze(SPEC(:,:,1))))), max(max(real(squeeze(SPEC(:,:,1))))) ]);
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    hold off;
    grid on
    pause(0.1)
end
%close(v)
end

if flag.plotsignals == 1
figure(3)
for i = 1:N3
    subplot(221)
    imagesc(-omega_4,omega_2,squeeze(real(I4WMt(:,:,i)))); hold on;
    contourf(-omega_4,omega_2,squeeze(real(I4WMt(:,:,i))),20);
    axis square
    title(sprintf('T = %4.2f ps', tau_3(i)/3e-2))
    xlabel('-\omega_t/cm^{-1}');
    ylabel('\omega_\tau/cm^{-1}');
    set(gca,'FontSize',20,'YDir','normal','Linewidth',1.5,'Box','off')
    %set(gcf,'color','w');
    colorbar
    title('Total')
    %colormap(brewermap([],'RdBu'))
    colormap(cmap)
    caxis( [min(min(real(squeeze(I4WMt(:,:,1))))), max(max(real(squeeze(I4WMt(:,:,1))))) ]);
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    hold off;
    grid on
    subplot(222)
    imagesc(-omega_4,omega_2,squeeze(real(ESE(:,:,i)))); hold on;
    contourf(-omega_4,omega_2,squeeze(real(ESE(:,:,i))),20);
    axis square
    title(sprintf('T = %4.2f ps', tau_3(i)/3e-2))
    xlabel('-\omega_t/cm^{-1}');
    ylabel('\omega_\tau/cm^{-1}');
    set(gca,'FontSize',20,'YDir','normal','Linewidth',1.5,'Box','off')
    %set(gcf,'color','w');
    colorbar
    title('ESE')
    %colormap(brewermap([],'RdBu'))
    colormap(cmap)
    caxis( [min(min(real(squeeze(ESE(:,:,1))))), max(max(real(squeeze(ESE(:,:,1))))) ]);
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    hold off;
    grid on
    subplot(223)
    imagesc(-omega_4,omega_2,squeeze(real(GSB(:,:,i)))); hold on;
    contourf(-omega_4,omega_2,squeeze(real(GSB(:,:,i))),20);
    axis square
    title(sprintf('T = %4.2f ps', tau_3(i)/3e-2))
    xlabel('-\omega_t/cm^{-1}');
    ylabel('\omega_\tau/cm^{-1}');
    set(gca,'FontSize',20,'YDir','normal','Linewidth',1.5,'Box','off')
    %set(gcf,'color','w');
    colorbar
    %colormap(brewermap([],'RdBu'))
    colormap(cmap)
    caxis( [min(min(real(squeeze(GSB(:,:,1))))), max(max(real(squeeze(GSB(:,:,1))))) ]);
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    title('GSB')
    hold off;
    grid on
    subplot(224)
    imagesc(-omega_4,omega_2,squeeze(real(ESA(:,:,i)))); hold on;
    contourf(-omega_4,omega_2,squeeze(real(ESA(:,:,i))),20);
    axis square
    title(sprintf('T = %4.2f ps', tau_3(i)/3e-2))
    xlabel('-\omega_t/cm^{-1}');
    ylabel('\omega_\tau/cm^{-1}');
    set(gca,'FontSize',20,'YDir','normal','Linewidth',1.5,'Box','off')
    %set(gcf,'color','w');
    colorbar
    %colormap(brewermap([],'RdBu'))
    colormap(cmap)
    caxis( [min(min(real(squeeze(ESA(:,:,1))))), max(max(real(squeeze(ESA(:,:,1))))) ]);
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    hold off;
    grid on
    title('ESA')
    pause(0.1)
end
%close(v)
end

if flag.scroll_spectrum == 1
    scrollImagesc(-omega_4,omega_2,real(SPEC));
end

if flag.rhomat == 1

clr = jet(Ne^2-Ne);

figure(4)

rhoM = zeros(Ne^2,N3);
rho_exc = zeros(Ne,Ne,N3);
rho_site = zeros(Ne,Ne,Ne);

rho0 = V_S'*rho0_site*V_S; %exciton basis
for i = 1:length(tau_3)
    rhoM(:,i) = expm(Red_exc*tau_3(i))*rho0(:);
    rho_exc(:,:,i) = reshape(rhoM(:,i),Ne,Ne);
    rho_site(:,:,i) = V_S*reshape(rhoM(:,i),Ne,Ne)*V_S';
end

         
legpop = cell(1,Ne);
legcoh = cell(1,Ne^2 - Ne);
legspec = cell(1,Ne);
c = 1;
for i = 1:Ne
    for j = 1:Ne
        ind = ind2D(Ne,i,j);
        % site basis (system)
        if i == j
            legpop{i} = [num2str(i),num2str(i)];
            subplot(2,2,1);
            plot(taus,abs(squeeze(rho_site(i,i,:))),'b.-','color',pltc(i,:)); hold on
            % plot(real(squeeze(rho_site(1,1,:)) + squeeze(rho_site(2,2,:))),'g')
            title('Populations: Site Basis')
            grid on
            axis square
            xlabel('\tau_2/fs')
            ylabel('\rho_{aa}(\tau_2)')
            set(gca,'FontSize',16,'Linewidth',2,'Box','off')

            subplot(2,2,3);
            plot(taus,abs(squeeze(rho_exc(i,i,:))),'b.-','color',pltc(i,:)); hold on
            % plot(real(squeeze(rho_site(1,1,:)) + squeeze(rho_site(2,2,:))),'g')
            title('Populations: Exciton Basis')
            grid on
            axis square
            xlabel('\tau_2/fs')
            ylabel('\rho_{aa}(\tau_2)')
            set(gca,'FontSize',16,'Linewidth',2,'Box','off')
        else
            legcoh{c} = [num2str(i),num2str(j)];
            subplot(2,2,2);
            plot(taus,abs(squeeze(rho_site(i,j,:))),'b.-','color',pltc(i,:)); hold on
            % plot(real(squeeze(rho_site(1,1,:)) + squeeze(rho_site(2,2,:))),'g')
            title('Coherences: Site Basis')
            grid on
            axis square
            xlabel('\tau_2/fs')
            ylabel('\rho_{ab}(\tau_2)')
            set(gca,'FontSize',16,'Linewidth',2,'Box','off')

            legcoh{c} = [num2str(i),num2str(j)];
            subplot(2,2,4);
            plot(taus,abs(squeeze(rho_exc(i,j,:))),'b.-','color',pltc(i,:)); hold on
            % plot(real(squeeze(rho_site(1,1,:)) + squeeze(rho_site(2,2,:))),'g')
            title('Coherences: Exciton Basis')
            grid on
            axis square
            xlabel('\tau_2/fs')
            ylabel('\rho_{ab}(\tau_2)')
            set(gca,'FontSize',16,'Linewidth',2,'Box','off')

            c = c+ 1;
        end
    end
end
subplot(2,2,1)
ll = legend(legpop);    set(ll,'FontSize',13,'Location','NorthEast');
subplot(2,2,3)
ll = legend(legpop);    set(ll,'FontSize',13,'Location','NorthEast');
subplot(2,2,2)
ll = legend(legcoh);    set(ll,'FontSize',13,'Location','NorthEast');
subplot(2,2,4)
ll = legend(legcoh);    set(ll,'FontSize',13,'Location','NorthEast');
hold off

figure(7)
hold on;
for i = 1:Ne
    plot(omega_bath, JJ{i}, 'Linewidth',2)
    title('Spectral Densities')
    grid on
    axis square
    xlabel('\omega_{bath}/cm^{-1}')
    ylabel('J_{ij}(\omega)')
    set(gca,'FontSize',16,'Linewidth',2,'Box','off')
    legspec{i} = [num2str(i), num2str(i)];
    pause(0.01)
end
ll = legend(legspec); set(ll,'FontSize',13,'Location','NorthEast');
hold off;

[~,inddec ] = min( abs( squeeze(rho_site(3,3,:)) - 1/exp(1) ) );
taudec = taus(inddec)

end

if flag.plot_time_trace == 1
    subplot(221)
    plot(tau_3,real(squeeze(mean(mean(SPEC,1),2))),'Linewidth',2)
    grid on
    xlabel('\tau_2')
    ylabel('Intensity')
    set(gca,'FontSize',17)
    subplot(222)
    plot(omega_2,abs(squeeze(mean(mean(SPEC,2),3))),'Linewidth',2)
    xlabel('\omega_2')
    ylabel('Intensity')
    set(gca,'FontSize',17)
    grid on
    subplot(223)
    plot(omega_4,abs(squeeze(mean(mean(SPEC,1),3))),'Linewidth',2)
    xlabel('\omega_4')
    ylabel('Intensity')
    set(gca,'FontSize',17)
    grid on
    subplot(224)
    plot(omega_3,abs(squeeze(mean(mean(fftshift(fft(SPEC,N3,3),3),1),2))),'Linewidth',2)
    xlabel('\omega_3')
    ylabel('Intensity')
    set(gca,'FontSize',17)
    grid on
end

%% Derived pump probe spectrum
Ipp = squeeze(mean(Stot,1));
hold on
for i = 1:50:N3
    plot(omega_4,Ipp(:,i));
end
hold off

%% Plot Intensities of Xpeaks
[~,E20] = min( abs( omega_2 - bandgap ) );
[~,E22] = min( abs( omega_2 - (bandgap+eegap) ) );
[~,E40] = min( abs( -omega_4 - bandgap) );
[~,E44] = min( abs( -omega_4 - (bandgap+eegap) ) );

Xloc = [E40, E44];
Yloc = [E20, E22];

count = 1;
for i = 1:numstates-1
    for j = 1:numstates-1
        subplot(2,2,count)
        plot(taus,squeeze(I4WMt(Yloc(i),Xloc(j),:)),'Linewidth',2)
        xlabel('\tau_2 [fs]')
        ylabel('Intensity')
        grid on
        set(gca,'Fontsize',15,'Linewidth',1.5,'Box','off')
        ll = legend(sprintf('(%4.0f,%4.0f)',omega_2(Yloc(i)) - bandgap, -omega_4(Xloc(j)) - bandgap ));
        set(ll,'FontSize',14,'Location','NorthEast');
        axis([taus(1),taus(end),0,inf])
        count = count + 1;
    end
end

%% Step 1: Background subtraction
mm = 1;
tauroi = tau_3(mm:end);
n3 = length(tauroi);
[U,s,v] = svd(reshape(real(Stot(:,:,mm:end)),N2*N4,n3),'econ');
K = length(find(diag(s) >= 1e-6));
lb = -100*ones(1,K);
ub = 10*ones(1,K);
x0 = (ub-lb).*rand(1,K) + lb;
% use global on R3dyn because coherences are much stronger in simulation
[xsolve,Gmatpop,Lmatpop,fitpop,fval,residual,jacobian] = myGlobalAnalysis_v2(U,s,v,tauroi, K, lb,ub, x0,'time',1,10000);
gammapop_fit = xsolve;
%% Examine population fit
close all
for i = 1:K
    figure(1)
    xx = Lmatpop(i,:);
    subplot(3,3,i)
    plot(tauroi,real(xx))
    grid on
    figure(2) 
    subplot(3,3,i)
    gm = real(reshape(Gmatpop(:,i),N2,N4));
    imagesc(omega_4,omega_2,gm); hold on;
    contourf(omega_4,omega_2,gm,20);
    plot(omega_4,-omega_4,'k--','Linewidth',2)
    title(sprintf('g=%4.2f',gammapop_fit(i)));
    xlabel(sprintf('%4.0f',i),'FontSize',16)
    set(gca,'FontSize',16,'Ydir','normal')
    axis square
    colorbar
    colormap(jet_white)
end
%% Subtract out baseline DC
close all
baselinepk = [1,4];
%GMbl = real(reshape(I4WMt,N2*N4,N3))*pinv(ones(1,N3));
%LL = GMbl*ones(1,N3);
%TML = exp(gammapop_fit(baselinepk)*tau_3);
%LL = Gmatpop(:,baselinepk)*TML;
Ydiff = real(Stot) - reshape(Gmatpop(:,baselinepk)*Lmatpop(baselinepk,:),N2,N4,N3);
plot(taus,real(squeeze(mean(mean(Ydiff,1),2))),'Linewidth',2)
xlabel('\tau_3 [fs]')
ylabel('Intensity')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
grid on
title('Baseline DC Subtracted')

%% Pick out good DC peaks and create DAS
poppeak = [1,4];
DAS = reshape(Gmatpop(:,poppeak)*Lmatpop(poppeak,:),N2,N4,N3);
plot(taus, squeeze(mean(mean(DAS,1),2) ) ,'Linewidth',2)
xlabel('\tau_3 [fs]')
ylabel('Intensity')
title('DAS fit Residual')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
grid on
%% Step 2: Global analysis on full spectrum fit to complex exponentials
mm = 1;
tauroi = tau_3(mm:end);
n3 = length(tauroi);
%ROI = Ydiff(:,:,mm:end);% - CAS(:,:,mm:end)
[U,s,V] = svd(reshape(Ydiff,N2*N4,n3),'econ');
%
K = length(find(diag(s) >= 1e-8));
[xsolve,Gmatcoh,Lmatcoh,fit,fval] = myGlobalAnalysis_v2(reshape(Ydiff,N2*N4,n3), tauroi, K, [], [], 0, 'time2', 0);
omega_fit = real(xsolve);
gamma_fit = imag(xsolve);
%% Examine fit - coherences
%omega_4roi  = omega_4;
%omega_2roi = omega_2;
close all
for i = 1:K
    figure(1)
    xx = Lmatcoh(i,:);
    subplot(5,5,i)
    plot(tauroi,real(xx))
    grid on
    figure(2)
    subplot(5,5,i)
    gm = real(reshape(Gmatcoh(:,i),N2,N4));
    imagesc(omega_4,omega_2,gm); hold on;
    contourf(omega_4,omega_2,gm);
    title(sprintf('w=%4.2f/g=%4.2f',omega_fit(i),gamma_fit(i)));
    xlabel(sprintf('%4.0f',i),'FontSize',16)
    axis square
    colorbar
    colormap(cmap)
end
%% Step 3: Pick out the good beat maps and make the coherence spectrum
close all
indcoh = [2,4];
dcpeak = [];
CAS = CAS + reshape(Gmatcoh(:,[indcoh,dcpeak])*Lmatcoh([indcoh,dcpeak],:),N2,N4,n3);
CASft = fftshift(fft(CAS,N3,3),3);
plot(omega_3, squeeze(mean(mean(CASft,1),2) ) ,'Linewidth',2)
xlabel('\omega_3 [cm^{-1}]')
ylabel('Intensity')
title('CAS')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off')
grid on
%% Plot the CAS and DAS
close all
for i = 1:N3
subplot(1,2,1)
contourf(-omega_4,omega_2,real(CAS(:,:,i)),20)
xlabel('-\omega_4 [cm^{-1}]')
ylabel('\omega_2 [cm^{-1}]')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off','Ydir','normal')
title('CAS')
axis square
grid on
colorbar
caxis( [min(min(real(squeeze(CAS(:,:,1))))), max(max(real(squeeze(CAS(:,:,1))))) ]);
colormap(jet_white)
ll = legend({num2str(taus(i)/1000) 'ps'}); set(ll,'FontSize',15,'Location','North');
subplot(1,2,2)
contourf(-omega_4,omega_2,real(DAS(:,:,i)),20)
xlabel('-\omega_4 [cm^{-1}]')
ylabel('\omega_2 [cm^{-1}]')
set(gca,'FontSize',16,'Linewidth',1.5,'Box','off','Ydir','normal')
title('DAS')
axis square
grid on
colorbar
caxis( [min(min(real(squeeze(DAS(:,:,1))))), max(max(real(squeeze(DAS(:,:,1))))) ]);
colormap(jet_white)
ll = legend({num2str(taus(i)/1000) 'ps'}); set(ll,'FontSize',15,'Location','North');
%hh = suptitle(sprintf('T = %4.0f fs',taus(i))); set(hh,'FontSize',16);
pause(0.1)
end

%% Run modefinder
specres = abs(omega_3(2) - omega_3(1));
cutoff = 200;
sel = (max(peakspec)-min(peakspec))/5000;
thresh = 0;
[omega_out,vout,n_modes0,fout] = omegaT_select3(omega_3,peakspec,sel,thresh,specres,cutoff);
omega_T = omega_out;

%omega_T = [0, 215, 430, 590, 645, 800];
%omega_T = [0,215,430];
%% Step 5: Dipole matrix reconstruction - single manifold basis spectra

% omega_T = [0,215,430,510];

% first to slow pathways (ZQC: 10 cm^-1)
%indpeak = [3,4,6,8,10];
indpeak = [zqcpeak,dcpeak];
omega_peak = real(xsolve([indpeak]));
Gamma_sqc_rect = Gamma_sqc_ge1;
lineshape_rect = 'lorentzian';

L = length(indpeak);
M = length(omega_T);

Qs = zeros(M^3,L);
Q_RS = 0;
chkbasis = 1;

% Reconstruction for R3
for l = 1:L
    
    X = Gmatcoh(:,indpeak(l)); % coherence beat map
  
    IM = zeros(N2,N4,M,M,M);   
   
    for a = 1:M
        for beta = 1:M
            for b = 1:M
    
                E_g0 = 0;
                E_ep_1 = abs(omega_T(a)) + omega_ge(3);
                E_gp = abs(omega_T(beta));
                E_ep_2 = abs(omega_T(b)) + omega_ge(3);
                E1 = E_ep_1 - E_g0;
                E2 = E_gp - E_ep_2;
                Eab = E_ep_1 - E_ep_2;
                % if the potential beating map entries actually correspond
                % to the coherence slice we are taking (this is
                % IMPORTANT!!!)
                % IMPORTANT: in spectrum, E3 = omega_a - omega_b and here
                % should also be omega_a - omega_b
                
                threshold =abs( Eab - omega_peak(l) );
                
                if threshold  < 2*specres              
                        I2 = myLineshape_v3(omega_2,E1,Gamma_sqc_rect,lineshape_rect);
                        I4 = myLineshape_v3(omega_4,E2,Gamma_sqc_rect,lineshape_rect);
                        IM(:,:,a,beta,b) = reshape(kron(I4,I2),[N2,N4]);
                end
            end
        end
    end
    
    IM_RS = reshape(IM,[N2*N4,M^3]);
    
    % taking real part loses no information because of Kramers-Kronig
    IM_RS = real(IM_RS);
    X = real(X);
    
    %figure(8)
    %subplot(2,2,l)
    %imagesc(omega_2,omega_4,real(squeeze(sum(sum(sum(IM,3),4),5)))); hold on;
    %contourf(omega_2,omega_4,real(squeeze(sum(sum(sum(IM,3),4),5)))); hold off;
    %colorbar
    %axis square
    %colormap(redblue_alt)

    
    %[U_basis,s_basis,V_basis] = csvd(IM_RS);
    % gcv works better
    %lambda = l_curve(U_basis,s_basis,X);
    %figure(2)
    %lambda = gcv(U_basis,s_basis,X); 
    %Q = tikhonov(U_basis,s_basis,V_basis,X,lambda);
    
    % regularization can introduce some weird noise for "perfect" fits
    Q = pinv(IM_RS)*X;
    Q_RS = Q_RS + reshape(Q,[M,M,M]);
    Qs(:,l) = Q;
    
    
    if chkbasis == 1
    xcheck = IM_RS*Q;   
    figure(4)
    subplot(221)
    imagesc(omega_2,omega_4,reshape(real(X),[N2,N4])); hold on;
    contourf(omega_2,omega_4,reshape(real(X),[N2,N4]),16); hold off;
    title('Beat Map'); xlabel('\omega_2'); ylabel('\omega_4');
    colormap(redblue_alt)
    colorbar
    axis('square')
    subplot(222)
    imagesc(omega_2,omega_4,real(reshape(xcheck,[N2,N4]))); hold on;
    contourf(omega_2,omega_4,real(reshape(xcheck,[N2,N4])),16); hold off;
    title('Reconstruction'); xlabel('\omega_2'); ylabel('\omega_4');
    colormap(redblue_alt)
    colorbar
    axis('square')
    subplot(223)
    plot(real(Q),'r.-')
    axis('square')
    title('Q')
    grid on
    pause(2)

    end
    
    
end

subplot(224)
plot((imag(Q_RS(:))),'r.-');  
hold on
plot(real(Q_RS(:)))
hold off
title('All Q')
grid on
axis('square')
ll = legend('Imag','Real'); set(ll,'FontSize',10,'Location','Best')




%% Recover Individual Transition Moments + Dpv Regression (scale factor works!)


%n_modes = 1;
%n_quanta = 10;
Bv = (1 + 1e-6)*ones(1,n_modes);

options = optimoptions('LSQNONLIN','Display','iter-detailed','Diagnostics','on');
options.FunctionTolerance = 1e-10;
options.OptimalityTolerance = 1e-12;
options.StepTolerance = 1e-14;

H = zeros(M^3);

for a = 1:M    
    for beta = 1:M
        for b = 1:M
            
            v1 = zeros(1,M); v1(1,a) = 1;
            v2 = zeros(M,M); v2(beta,a) = 1;
            v3 = zeros(M,M); v3(beta,b) = 1;
            v4 = zeros(1,M); v4(1,b) = 1;
            
            Ytemp = kron(v4',(kron(v3,kron(v2',v1))));
            [I,J] = find(Ytemp == 1);
            H(I,J) = Q_RS(a,beta,b);
            
            %if a == b
            %    H(I,J) = 0;
            %end
            
        end
    end
end

x0 = rand(M,M); lb = min(diag(H))*ones(M); ub = max(diag(H))*ones(M);
fun = @(x)diag(kron(x(1,:)',(kron(x,kron(x',x(1,:)))))) - diag(H);
[xsolve1,fval1] = lsqnonlin(fun,x0,lb,ub,options);
xsolve1

lb = [zeros(1,n_modes),1e-9];
ub = [ones(1,n_modes),1];
x0 = (ub-lb).*rand(1,n_modes+1) + lb;

fun = @(x) myFCRegression_v4(Bv, x, M, abs(xsolve1), n_modes, n_quanta, []);

[xsolve2,fval2] = lsqnonlin(fun,x0,lb,ub,options);
xsolve2
fval2


%% Step 6: Minimize error to analytical FC matrix formula to extract displacements (and curvatures)

% for R3 + R4, Qs -> Qs/2 to make this work since dipole moment is same for
% both pathways
options = optimset('Display','Iter-Detailed','TolFun',1e-8,'TolX',1e-8);
x0 = rand(1,n_modes);
fun = @(x) myFCRegression(Bv, x, M, Qs, n_modes, n_quanta,[]);
[xsolve2,fval2] = fminsearch(fun,x0,options); 
xsolve2
Dpv
