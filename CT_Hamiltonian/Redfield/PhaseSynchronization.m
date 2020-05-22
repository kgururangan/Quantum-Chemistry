clear all
clc
close all
%
addpath('/Users/karthik/Documents/MATLAB/2DPES/QPT');
addpath('/Users/karthik/Documents/MATLAB/Colormaps');
addpath('/Users/karthik/Documents/MATLAB/Additional_Functions');
addpath('/Users/karthik/Documents/MATLAB/Additional_Functions/regu');

pltc = default_cmap();
%

secular = 0;

omega_c = 300;
gamma_i = 0.5;
J0 = 0.5*omega_c;
Delta = 2*omega_c;
theta = omega_c;

n_quanta = 1;
omega_e = [Delta, Delta];
omega_ge = [-Delta, Delta];
Dpv = [0.1, 0.1];

Cre = diag(sqrt([1:n_quanta-1]),-1);
Ann = Cre';
Num = Cre*Ann;
H_e  = cell(1,2);
for i = 1:2
    H_e{i} = omega_e(i)*Num + Dpv(i)*(Cre+Ann)*omega_e(i) + omega_ge(i)*eye(n_quanta);
end

H_S = [H_e{1}, J0*eye(n_quanta); J0*eye(n_quanta), H_e{2}];
M = length(H_S);

Nw = 512;
domega = 10; % cm^-1
omega = domega:domega:(Nw)*domega; % cm^-1
domega = abs(omega(2) - omega(1));

JJ{1,1} = gamma_i*omega.^3.*exp(-omega/omega_c)./(pi*omega_c.^2); % super-ohmic
JJ{2,2} = JJ{1,1};
JJ{1,2} = 0.1*JJ{1,1}; 
%JJ{1,2} = zeros(1,Nw);
JJ{2,1} = JJ{1,2};
[V_S,E_S] = eig(H_S); % E_S = V_S'*H_S*V_S
E = diag(E_S);

%
tic
c = 3e-5; % cm/fs
Nt = 2024;
dt_fs = 0.5; % fs
dt = 2*pi*dt_fs*c; % rad*cm

% t and tau are same axis
t = 0:dt:(Nt-1)*dt; % rad*cm
t = t';
%
taus = t*omega_c; 
%
n = @(w) (1 + exp(omega/theta) )^-1;

Ym = zeros(M,M,Nw,Nt);

for i = 1:M
    for j = 1:M
        for k = 1:Nw
            
            %Ym(i,j,k,:) = deltafnc(i,j)*( JJ(k)*(cos(omega(k)*t)*coth(omega(k)/2/theta) - 1i*sin(omega(k)*t)) );
            Ym(i,j,k,:) =  JJ{i,j}(k)*(cos(omega(k)*t)*coth(omega(k)/2/theta) - 1i*sin(omega(k)*t)) ;
        end
    end
end

Qm = zeros(M,M,M,M,Nt);

for alpha = 1:M
    for beta = 1:M
        for gamma = 1:M
            for delta = 1:M 
                % sum over states' spectral densities J_ij
                temp1 = 0;
                for n = 1:M
                    for m = 1:M
                        
                       temp1  = temp1 + conj(V_S(n,alpha))*V_S(n,beta)*...
                            conj(V_S(m,gamma))*V_S(m,delta)*squeeze(sum(Ym(n,m,:,:),3))*domega;                                    
                    end
                end      
                Qm(alpha,beta,gamma,delta,:) = temp1;
            end
        end
    end
end

Rtot = zeros(M,M,M,M);
Rtot_NS = cell(M,M,M,M);

for alpha = 1:M
    for beta = 1:M
        for gamma = 1:M
            for delta = 1:M
                
                R1temp = 0;
                for eta = 1:M
                    R1temp = R1temp + sum(exp(-1i*(E(eta) - E(gamma))*t).*...
                        squeeze(Qm(gamma,eta,eta,alpha,:))); 
                end
                R1v = deltafnc(beta,delta)*R1temp;  
                
                R2v = sum(exp(-1i*(E(delta) - E(beta))*t).*...
                    squeeze(conj(Qm(delta,beta,gamma,alpha,:))));
                
                R3v = sum(exp(-1i*(E(alpha) - E(gamma))*t).*...
                    squeeze(Qm(gamma,alpha,delta,beta,:)));

                R4temp = 0;
                for eta = 1:M
                    R4temp = R4temp + sum(exp(-1i*(E(delta) - E(eta))*t).*...
                        squeeze(conj(Qm(delta,eta,eta,beta,:))));
                end
                R4v = deltafnc(gamma,alpha)*R4temp; 
                
                Rtot(alpha,beta,gamma,delta) = -1i*(E(gamma) - E(delta))*...
                    deltafnc(alpha,gamma)*deltafnc(delta,beta) - (R1v - R2v - R3v + R4v)*dt;  
                
                %Rtot_NS{alpha,beta,gamma,delta} = @(time) -1i*(E(gamma) - E(delta))*...
                %    deltafnc(alpha,gamma)*deltafnc(delta,beta) - exp(1i*(alpha - beta - gamma + delta)*time)*(R1v - R2v - R3v + R4v)*dt;
                
                 Rtot_NS{alpha,beta,gamma,delta} = @(time) -1i*(E(gamma) - E(delta))*...
                    deltafnc(alpha,gamma)*deltafnc(delta,beta) - exp(1i*(alpha - beta - gamma + delta)*time)*(R1v - R2v - R3v + R4v)*dt;
            end
        end
    end
end
toc
R_rs = reshape(Rtot,M^2,M^2)
%% Evaluate rho and X_i
rho0_site = zeros(M,M);
%rho0_site = [0 0 0 0;
%             0 1 0 0;
%             0 0 0 0;
%             0 0 0 0];
rho0_site = [0, 0; 0, 1];
rho0 = V_S'*rho0_site*V_S; %exciton basis
R_rs = reshape(Rtot,M^2,M^2); % relaxation superoperator (d(rh0)/dt = R_rs rho, rho = M^2 x Nt
[H,B] = eig(R_rs);
Hi = inv(H);

% Calculate rho and Xav
Bi = kron(Cre+Ann,eye(2)) + kron(eye(2),Cre+Ann);
Xav = zeros(M^2,Nt);

% Evaluate in Time
rho_site = zeros(M,M,Nt);
rhoM = zeros(M^2,Nt);

W = zeros(4,Nt,2);

if secular == 1
    for i = 1:length(t)
        rhoM(:,i) = expm(R_rs*t(i))*rho0(:);
        rho_site(:,:,i) = V_S*reshape(rhoM(:,i),M,M)*V_S'; 
        Temp = rho_site(:,:,i)*Bi;
        Xav(:,i) = Temp(:);
    end
else % non-secular RK4 integration
     
    dteval = dt; % rad*cm
    Nteval = Nt;
    teval = 0:dteval:(Nteval-1)*dteval; % rad*cm
    Kt = zeros(M^2,M^2,Nteval);
    
    for k = 1:Nteval
        Temp = zeros(M,M,M,M);
        for alpha = 1:M
            for beta = 1:M
                for gamma = 1:M
                    for delta = 1:M
                        Temp(alpha,beta,gamma,delta) = Rtot_NS{alpha,beta,gamma,delta}(teval(k));
                    end
                end
            end
        end
        Kt(:,:,k) = reshape(Temp,M^2,M^2);
        W(:,k,1) = imag(Kt(2,:,k)); W(:,k,2) = real(Kt(2,:,k));
    end
    
    %Kt_rs = reshape(Kt,Nteval,M^4);
    Kt_rs = permute(Kt,[3,1,2]);
    
    L = @(ti) interp1(teval,Kt_rs,ti);
    
    x0 = rho0(:);
    rhoM = zeros(M^2,Nteval);
    rho_site = zeros(M,M,Nteval);
    
    for i = 1:Nteval
        
        m1 = squeeze(L(teval(i)));
        m2 = squeeze(L(teval(i)+0.5*dteval));
        m3 = squeeze(L(teval(i)+0.5*dteval));
        m4 = squeeze(L(teval(i)+dteval));
        
        % Fourth-order Runge-Kutta (RK4)
        k1 = dteval*m1*x0;
        k2 = dteval*m2*(x0+0.5*k1);
        k3 = dteval*m3*(x0+0.5*k2);
        k4 = dteval*m4*(x0+k3);
        
        x0 = x0 + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        
        rhoM(:,i) = x0;
        rho_site(:,:,i) = V_S*reshape(rhoM(:,i),M,M)*V_S';

    end
    
end


for i = 1:4
hold on
subplot(121)
plot(t,squeeze(W(i,:,1)),'Linewidth',2)
hold off
axis square
grid on
hold on
subplot(122)
plot(t,squeeze(W(i,:,2)),'Linewidth',2)
hold off
axis square
grid on
end

%% Plotting 
close all
legpop = cell(1,M);
legcoh = cell(1,M^2-M);
c = 1;
for i = 1:M
    for j = 1:M
        ind = ind2D(M,i,j);
        % site basis (system)
        if i == j
            legpop{i} = [num2str(i),num2str(i)];
            subplot(2,2,1);
            plot(taus,abs(squeeze(rho_site(i,i,:))),'b.-','color',pltc(i,:)); hold on
            % plot(real(squeeze(rho_site(1,1,:)) + squeeze(rho_site(2,2,:))),'g')
            title('Populations')
            grid on
            axis square
            axis tight
            xlabel('\tau_2 [1/\omega_c]')
            ylabel('\rho_{pop}(\tau_2)')
            set(gca,'FontSize',16,'Linewidth',2,'Box','off')
        else
            legcoh{c} = [num2str(i),num2str(j)];
            subplot(2,2,2);
            plot(taus,abs(squeeze(rho_site(i,j,:))),'b.-','color',pltc(i,:)); hold on
            % plot(real(squeeze(rho_site(1,1,:)) + squeeze(rho_site(2,2,:))),'g')
            title('Coherences')
            grid on
            axis square
            axis tight
            xlabel('\tau_2 [1/\omega_c]')
            ylabel('\rho_{pop}(\tau_2)')
            set(gca,'FontSize',16,'Linewidth',2,'Box','off')
            c = c+ 1;
        end
    end
end
subplot(2,2,1)
ll = legend(legpop);    set(ll,'FontSize',13,'Location','NorthEast');
subplot(2,2,2)
ll = legend(legcoh);    set(ll,'FontSize',13,'Location','NorthEast');
hold off


VSS = kron(V_S,V_S);
VSST = kron(V_S', V_S');

posv = zeros(M,M);
for a = 1:M
    for c = 1:M
        
        b = a; d = c;
        I1 = ind2D(M,a,b);
        I2 = ind2D(M,c,d);
        
        posv(c,a) = ind2D(M^2,I1,I2);
    end
end
pos = posv(:); % vector containing population terms for 1 x Ne^4 arrangement

subplot(224)
TMred = TM_matrix_v2(VSS*R_rs*VSST,t,0,'expm');
plot(taus,TMred(pos,:),'Linewidth',2)
grid on
axis square
axis tight
xlabel('\tau_2 [1/\omega_c]')
ylabel('T_{abcd}^{site}(\tau_2)')
title('Process Matrix')
set(gca,'FontSize',16,'Linewidth',2,'Box','off')
    
%% Oscillator Expectation values
figure(10)
hold on
for i = 1:M
    J = ind2D(M,i,i);
    plot(taus,Xav(J,:),'Linewidth',2)
    legpop{i} = [num2str(i),num2str(i)];
end
grid on
xlabel('\tau_2 [1/\omega_c]')
ylabel('\langle X _i\rangle(\tau_2)')
set(gca,'FontSize',16,'Linewidth',2,'Box','off')

i1 = 3;
i2 = 4;
x1 = real(Xav(i1,:));
x2 = real(Xav(i2,:));
di = 1; dti = di/(2*pi*3e-5);
C = zeros(1,Nt-di);
xC = zeros(1,length(C));
for i = 1:Nt-di
    f1 = x1 - sum(x1(i:i+di))*dti;
    f2 = x2 - sum(x2(i:i+di))*dti;
    
    %[co] = corrcoef(f1,f2);
    %C(i) = co(1,2);
    %xC(i) = mean([i:i+di]*dt/(2*pi*3e-5));
end
plot(xC,C,'Linewidth',2)
grid on
xlabel('\tau_2 / fs')
ylabel('C(t)')
set(gca,'FontSize',16,'Linewidth',2,'Box','off')

legpop{M+1} = 'C(x_1,x_2)';
ll = legend(legpop);    set(ll,'FontSize',13,'Location','NorthEast');

hold off
%%
x1 = sin(t*200);
x2 = sin(t*200+2*pi);
r = zeros(1,Nt);
di = 1;
C = zeros(1,Nt-di);
lags = [1:Nt-di];
for i = 1:Nt-di
    f1 = (x1(i:i+di)-1/(di*dt)*trapz(x1(i:di)));
    f2 = (x2(i:i+di)-1/(di*dt)*trapz(x2(i:i+di)));
    Temp = trapz(f1.*f2);
    Temp2 = sqrt(trapz(f1.^2.*f2.^2));
    C(i) = Temp/Temp2;
end
plot(t,x1,t,x2,t,r)

%%
close all
clc
clear all

load('nir_shootout_2002.mat')

% training data

% 1 and 2 are same dataset but taken with different spectrometer
lambda1 = calibrate_1.axisscale{2,1};
IMtrain1 = calibrate_1.data;

lambda2 = calibrate_2.axisscale{2,1};
IMtrain2 = calibrate_2.data;

ydescr = calibrate_Y.axisscale{2,1};
yTrain = calibrate_Y.data; % (:,1) = weight, (:,2) = hardness, (:,3) = assay

% testing data

lambda3 = test_1.axisscale{2,1};
IMtest1 = test_1.data;

lambda4 = test_2.axisscale{2,1};
IMtest2 = test_2.data;

yTest = test_Y.data;

% validation data

lambda5 = validate_1.axisscale{2,1};
IMval1 = validate_1.data;

lambda6 = validate_2.axisscale{2,1};
IMval2 = validate_2.data;

yVal = validate_Y.data;

% lambda 1-6 are all the same wavelength axis

%save('nir_shootout_2002_kg.mat')


    