function [Red_exc, GammaPop, JJ, Q] = RedfieldTensor(x, domega, omega_max, npts_per, tau_max, E, V_S, kbT)

% I don't htink we need to add GammaM, this is relaxation > dephasing by
% definition


theta = kbT;
% N3 = length(tau_3);
M = length(E); % number of excited state for ESE
% tau_3 = reshape(tau_3,N3,1);
% t = tau_3;
% Nt = length(tau_3);
% dt = abs(tau_3(2)-tau_3(1));
%Ntr = length(Itr); % number of states transferring via spectral density


% solvent bath omega coordinate
%domega = 10;
%omega_max = 5000;
%omega = 10:domega:omega_max;
omega = domega:domega:omega_max;
Nw = length(omega);
%omega(1) = 1e-6;

omega_c = zeros(M,1);
gamma_i = zeros(M,1);

JJ = cell(M,Nw);
for i  = 1:M
        s = 2*(i-1);
        omega_c(i) = x(s+1);
        gamma_i(i) = x(s+2);
        if omega_c(i) ~= 0
            JJ{i} =  gamma_i(i)*omega.^3.*exp(-omega/omega_c(i))./(pi*omega_c(i).^2); % super-ohmic
        else
            JJ{i} = zeros(1,Nw);
        end
end


tau_unit = 1/omega_c(1);
dt = 1/npts_per*tau_unit;
t = 0:dt:tau_max*tau_unit;
t = t';

Nt = length(t);
 
% plot(omega,JJ{1})

%omega_c = 1000;
%gamma_i = 0.5;
%Gammadab = 0;

%c = 3e-5; % cm/fs

%omega_3 = reshape(omega_3,N3,1) + max(omega_3); % doesn't work!
%omega_3 = ifftshift(omega_3);

%domega = abs(omega_3(2)-omega_3(1));

%omega = linspace(0,2000,Nw);

%omega = omega_3;

%domega = abs(omega(2) - omega(1));

%Ha = diag([2*omega_c,-2*omega_c]);

%GammaM = zeros(M^2,M^2);

%Jp = gamma_i*omega.^3.*exp(-omega/omega_c)./(pi*omega_c.^2); % super-ohmic
%Jp = (omega >= 0).*Jp;
%Jn = gamma_i*(-omega_3).^3.*exp(omega_3/omega_c)./(pi*omega_c.^2); % super-ohmic
%Jn = (omega_3 < 0).*Jn;
%J = Jn + Jp;
%J = Jp + flipud(Jp);
%J = Jp;

%J0 = 0.5*omega_c;

% states
%Jmat = zeros(size(Ha));
%Jmat(Itr,Itr) = J0*(ones(Ntr)-eye(Ntr));

%H_S = Ha + Jmat;
%H_S
%H_S= Ha( + J0*(ones(M) - eye(M));
%H_S = [2*omega_c, J0; J0, -2*omega_c];

%[V_S,E_S] = eig(H_S); % E_S = V_S'*H_S*V_S
%E = diag(E_S);


% bath correlation function <F(t)F(0)>_b = <(-kx(t))(-kx(0))> = spring
% force
Y = zeros(M,M,Nw,Nt);

for i = 1:M
    for j = 1:M
        for k = 1:Nw
            
            % harmonic bath model used here
            % deltafcn(i,j) = no cross-correlations for site fluctuations on
            % i,j coupled to bath
            J = JJ{i};
            Y(i,j,k,:) = deltafnc(i,j)*( J(k)*(cos(omega(k)*t)*coth(omega(k)/2/theta) - 1i*sin(omega(k)*t)) );
            %Y(i,j,k,:) = ( J(k)*(cos(omega(k)*t)*coth(omega(k)/2/theta) - 1i*sin(omega(k)*t)) );

        end
    end
end


Q = zeros(M,M,M,M,Nt);

for alpha = 1:M
    for beta = 1:M
        for gamma = 1:M
            for delta = 1:M
                
                % sum over states' spectral densities J_ij
                temp = 0;
                for n = 1:M
                    for m = 1:M
                        
                       temp  = temp + conj(V_S(n,alpha))*V_S(n,beta)*...
                            conj(V_S(m,gamma))*V_S(m,delta)*squeeze(sum(Y(n,m,:,:),3))*domega;
                       %temp  = temp + conj(V_S(n,alpha))*V_S(n,beta)*...
                       %     conj(V_S(m,gamma))*V_S(m,delta)*trapz(omega,squeeze(Y(n,m,:,:)),1);
                                      
                    end
                end
                
                Q(alpha,beta,gamma,delta,:) = temp;

            end
        end
    end
end


Rtot = zeros(M,M,M,M);

for alpha = 1:M
    for beta = 1:M
        for gamma = 1:M
            for delta = 1:M
                
                R1temp = 0;
                for eta = 1:M
                    R1temp = R1temp + sum(exp(-1i*(E(eta) - E(gamma))*t).*...
                        squeeze(Q(gamma,eta,eta,alpha,:))); % sum(.) is integrating over tau
%                         squeeze(Q(alpha,eta,eta,gamma,:)));
                end
                R1v = deltafnc(beta,delta)*R1temp; 
                
                R2v = sum(exp(-1i*(E(delta) - E(beta))*t).*...
                    squeeze(conj(Q(delta,beta,gamma,alpha,:))));
%                 squeeze(Q(alpha,gamma,delta,beta,:)));
                
                R3v = sum(exp(-1i*(E(alpha) - E(gamma))*t).*...
                    squeeze(Q(gamma,alpha,delta,beta,:)));
%                 squeeze(Q(delta,beta,alpha,gamma,:)));
                R4temp = 0;
                for eta = 1:M
                    R4temp = R4temp + sum(exp(-1i*(E(delta) - E(eta))*t).*...
                        squeeze(conj(Q(delta,eta,eta,beta,:))));
%                     squeeze(Q(eta,beta,delta,eta,:)));
                end
                R4v = deltafnc(gamma,alpha)*R4temp;
                
                Rtot(alpha,beta,gamma,delta) = -1i*(E(gamma) - E(delta))*...
                    deltafnc(alpha,gamma)*deltafnc(delta,beta) - (R1v - R2v - R3v + R4v)*dt;
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

%Red_exc = reshape(Rtot,M^2,M^2);
%Red_site = kron(V_S,V_S)*Red_exc*kron(V_S',V_S');

%Red_site_rs = kron(V_S,V_S)*Red_exc_rs*kron(V_S',V_S');

%Red_rs = reshape(Red,M^2,M^2);
%TMred = TM_matrix_v2(Red_rs,tau_3,0);
%plot(tau_3,real(TMred))

end

function xout = deltafnc(a,b)
    if a == b
        xout = 1;
    else
        xout = 0;
    end
end

