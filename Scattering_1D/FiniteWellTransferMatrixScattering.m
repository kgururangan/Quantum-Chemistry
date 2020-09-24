clear all
clc
close all

V0 = 8e-04;
a = 3;
m = 1822;
E = linspace(1.0e-05,1.0e-02,1000);

T = zeros(1,length(E));
R = zeros(1,length(E));

for i = 1:length(E)
    
    k = sqrt(2*m*E(i));
    kappa = sqrt(2*m*(E(i) - V0));

%     M1 = [1, 1;
%           k, -k];
%       
%     M2 = [1, 1;
%           kappa, -kappa];
%       
%     M3 = [exp(1i*kappa*a), exp(-1i*kappa*a);
%           kappa*exp(1i*kappa*a), -kappa*exp(-1i*kappa*a)];
%     
%     Mtilde = M3*(M2\M1);
%     
%     M = [Mtilde(1,2), -exp(1i*k*a);
%          Mtilde(2,2), -k*exp(1i*k*a)];
%      
%     b = [-Mtilde(1,1), -Mtilde(2,1)]';
% 
%     
%     x = M\b;
%     
%     R(i) = sqrt(conj(x(1))*x(1));
%     %T(i) = sqrt(1 - R(i)^2);
%     T(i) = sqrt(conj(x(2))*x(2));
    
    
    A = [1, -1, -1, 0;
         -k, -kappa, kappa, 0
         0, exp(1i*kappa*a), exp(-1i*kappa*a), -exp(1i*k*a); 
         0, kappa*exp(1i*kappa*a), -kappa*exp(-1i*kappa*a), -k*exp(1i*k*a)];
     
    b = [-1, -k, 0, 0]';

%     A = [exp(-1i*k*a), -exp(-1i*kappa*a), exp(1i*kappa*a), 0;
%         -k*exp(1i*k*a), -kappa*exp(-1i*kappa*a), kappa*exp(1i*kappa*a), 0;
%         0, exp(1i*kappa*a), exp(-1i*kappa*a), -exp(1i*k*a);
%         0, kappa*exp(1i*kappa*a), -kappa*exp(-1i*kappa*a), -k*exp(1i*k*a)];
%         
%     b = [-exp(-1i*k*a), -k*exp(-1i*k*a), 0, 0]';
%     
    x = A\b;
    
    R(i) = sqrt(conj(x(1))*x(1));
    T(i) = sqrt(conj(x(4))*x(4));



%     fprintf('E = %4.10f\n',E(i))
%     fprintf('|T|^2 = %4.10f\n',T(i)^2)
%     fprintf('|R|^2 = %4.10f\n\n',R(i)^2)
    
end

plot(E/V0,T*100,E/V0,R*100,'Linewidth',2)
xlabel('E / V_0')
ylabel('[%]')
set(gca,'FontSize',15,'Linewidth',2,'Box','off')
axis([-Inf, Inf, 0, 120])
grid on
ll = legend('Transmission','Reflectance'); set(ll,'FontSize',14,'Location','Best')

%% Transfer matrix method

clear all
clc
close all

V = 1.0;
V0 = 0.0;

a = 1;
b = 1;
m = 100;
E = linspace(1.5,3.5,1000);

T = zeros(1,length(E)); 
R = zeros(1,length(E)); 

%phi_left = [0,1]';

L = 20;

for i = 1:length(E)
    
    P = @(v,a) PropatationMatrix(v,a,E(i),m);
    D = @(v1,v2) DiscontinuityMatrix(v1,v2,E(i),m);
    
%     M = D(0,V0)*P(V0,a)*D(V0,0);
%     phi_right = M*phi_left;
%     T(i) = sqrt(conj(1/M(2,2))*(1/M(2,2)));
%     R(i) = sqrt(conj(M(1,2))*M(1,2));

    M = eye(2);
    V0 = V0; V1 = V;
    M1 = P(V0,b)*D(V0,V1)*P(V1,a)*D(V1,V0);
    
    M = M1^L;
%     for j = 1:L
%         M1 = P(V0,b)*D(V0,V1)*P(V1,a)*D(V1,V0);
%         M = M * M1;
%     end
    %phi_right = M*phi_left;
    T(i) = (conj(1/M(2,2))*(1/M(2,2)));
    R(i) = (conj(M(2,1))*M(2,1)/(conj(M(1,1))*M(1,1)));

%     fprintf('E = %4.10f\n',E(i))
%     fprintf('|T|^2 = %4.10f\n',T(i)^2)
%     fprintf('|R|^2 = %4.10f\n\n',R(i)^2)
    
end

plot(E/V1,T*100,'Linewidth',2)
xlabel('E / V_0')
ylabel('Transmission [%]')
set(gca,'FontSize',15,'Linewidth',2,'Box','off')
%axis([-Inf, Inf, 0, 120])
grid on



function P = PropatationMatrix(v,a,e,m)
    k = sqrt(2*m*(e - v));
    P = [exp(1i*k*a), 0;
            0       exp(-1i*k*a)];
end

function [D] = DiscontinuityMatrix(v1,v2,e,m)
    k1 = sqrt(2*m*(e-v1));
    k2 = sqrt(2*m*(e-v2));
    da = 0.5*(1 + k1/k2);
    db = 0.5*(1 - k1/k2);
    D = [da, db;
         db, da];
end