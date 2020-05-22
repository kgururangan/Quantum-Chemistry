clear all
%clc
close all

% E0 = 12500;
% Delta = 2000;
% J0 = 1000;
% 
% H = [E0, J0;
%      J0, E0+Delta];
%  

Delta = 1000;
omega_c = 2*Delta;
tau_unit = 1/omega_c;
J0 = 0.5*omega_c;
theta = omega_c;
gamma_i = 1;

npts_per = 60;
tau_max = 40;
dt = 1/npts_per*tau_unit;
t = 0:dt:tau_max*tau_unit;

t = reshape(t,length(t),1);

H = [-Delta, J0;
      J0, Delta];
  
[V_S, E] = eig(H);

M = size(H,1);


domega = 10; omega_max = 20*omega_c;
omega = domega:domega:omega_max;
%

Nt = length(t); Nw = length(omega);

J =  gamma_i*omega.^3.*exp(-omega/omega_c)./(pi*omega_c.^2); 

subplot(221)
plot(omega,J)

Y = zeros(M,M,Nw,Nt);
for i = 1:M
    for j = 1:M
        for k = 1:Nw
            Y(i,j,k,:) = deltafnc(i,j)*( J(k)*(cos(omega(k)*t)*coth(omega(k)/2/theta) - 1i*sin(omega(k)*t)) );
        end
    end
end

Q = zeros(M,M,M,M,Nt);

subplot(222)
hold on
for alpha = 1:M
    for beta = 1:M
        for gamma = 1:M
            for delta = 1:M
                
                temp = 0;
                for n = 1:M
                    for m = 1:M      
                       temp  = temp + conj(V_S(n,alpha))*V_S(n,beta)*...
                            conj(V_S(m,gamma))*V_S(m,delta)*squeeze(sum(Y(n,m,:,:),3))*domega;
                                      
                    end
                end
                
                Q(alpha,beta,gamma,delta,:) = temp;
                plot(t,squeeze(Q(alpha,beta,gamma,delta,:)))
            end
        end
    end
end
hold off

%
Rtot = zeros(M,M,M,M);

for alpha = 1:M
    for beta = 1:M
        for gamma = 1:M
            for delta = 1:M
                
                R1temp = 0;
                for eta = 1:M
                    R1temp = R1temp + sum(exp(-1i*(E(eta) - E(gamma))*t).*squeeze(Q(gamma,eta,eta,alpha,:)))*dt;
                end
                R1v = deltafnc(beta,delta)*R1temp; 
                
                R2v = sum(exp(-1i*(E(delta) - E(beta))*t).*squeeze(conj(Q(delta,beta,gamma,alpha,:))))*dt;
                
                R3v = sum(exp(-1i*(E(alpha) - E(gamma))*t).*squeeze(Q(gamma,alpha,delta,beta,:)))*dt;
                R4temp = 0;
                for eta = 1:M
                    R4temp = R4temp + ...
                        sum(exp(-1i*(E(delta) - E(eta))*t).*squeeze(conj(Q(delta,eta,eta,beta,:))))*dt;
                end
                R4v = deltafnc(gamma,alpha)*R4temp;
                
                Rtot(alpha,beta,gamma,delta) = -1i*(E(gamma) - E(delta))*...
                    deltafnc(alpha,gamma)*deltafnc(delta,beta) - (R1v - R2v - R3v + R4v);
                %Rtot(alpha,beta,gamma,delta) = (R1v - R2v - R3v + R4v)*dt;
                
                
            end
        end
    end
end


Red_exc = Rtot;

R = reshape(Red_exc,M^2,M^2);
GammaPop = zeros(M);
for p = 1:M
    for q = 1:M
        I1 = ind2D(M,p,p);
        I2 = ind2D(M,q,q);
        GammaPop(p,q) = R(I1,I2);
    end
end

GammaPop = real(GammaPop)

rho_site = zeros(M,M,length(t));
rho_exc = zeros(M,M,length(t));

rho_site(:,:,1) = [0,0;
                  0,1];
rho_exc(:,:,1) = V_S*squeeze(rho_site(:,:,1))*V_S';
for i = 2:length(t)
    rho_exc(:,:,i) = reshape(expm(R*t(i))*reshape(rho_exc(:,:,1),M^2,1),M,M);
    rho_site(:,:,i) = V_S*squeeze(rho_exc(:,:,i))*V_S';
end

for i = 1:M
    subplot(223)
    plot(t,squeeze(rho_exc(i,i,:)))
    hold on
    xlabel('t / cm')
end
hold off

for i = 1:M
    subplot(224)
    plot(t,squeeze(rho_site(i,i,:)))
    hold on
    xlabel('t / cm')
end
hold off

