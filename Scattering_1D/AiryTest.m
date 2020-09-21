%% Trapezoidal barrier (1 line)

clear all
clc
close all

m = 100;
alpha = 0.1;
a = 1;
V0 = 1;
V1 = 2;

E = linspace(1,5,1000) * V0;

T = zeros(1,length(E));
T_exact = zeros(1,length(E));
for i = 1:length(E)
    
    k = sqrt(2*m*E(i));
    k2 = conj(k)*k;
    
    l01 = (2*m*alpha)^(-1/3);
    z1 = @(x) 1/l01*(x - (E(i) - V0)/alpha);
    
    ML = [1, 1;
          1i*k, -1i*k];
      
    M1 = [airy(0,z1(0)), airy(2,z1(0));
         1/l01*airy(1,z1(0)), 1/l01*airy(3,z1(0))];
     
    M2 = [airy(0,z1(a)), airy(2,z1(a));
         1/l01*airy(1,z1(a)), 1/l01*airy(3,z1(a))];
     
    MR = [exp(1i*k*a), 0;
         1i*k*exp(1i*k*a), 0];
     
    %Mtot = inv(ML) * M1 * inv(M2) * MR;
    Mtot = (ML\M1)*(M2\MR);
   
    
    
    T(i) = 1/Mtot(1,1);
    T(i) = conj(T(i))*T(i);
    
    % exact result
    f = airy(0,z1(0));
    fp = 1/l01*airy(1,z1(0));
    g = airy(2,z1(0));
    gp = 1/l01*airy(3,z1(0));
    q = airy(0,z1(a));
    qp = 1/l01*airy(1,z1(a));
    r = airy(2,z1(a));
    rp = 1/l01*airy(3,z1(a));
    
    T_exact(i) = (2i*k*(r*qp - q*rp))/...
        ( (-k2*f*r - fp*rp + k2*g*q + gp*qp) + 1i*k*(fp*r - f*rp - gp*q + g*qp) );
    T_exact(i) = conj(T_exact(i))*T_exact(i);
    
end

yyaxis left
plot(E/V0,T_exact)
yyaxis right
plot(E/V0,T)


%% Triangular barrier (2 lines)

clear all
clc
close all

m = 100;
a = 1;
b = 1;

V0 = 1;
V1 = 1.1;

alpha = [(V1 - V0)/a, (V0 - V1)/b];

E = linspace(0,8,1000);

T = zeros(1,length(E));
for i = 1:length(E)
    
    k = sqrt(2*m*E(i));
    k2 = conj(k)*k;
    
    l0 = airylengthscale(alpha,m);
    z1 = @(x) 1/l0(1)*(x - (E(i) - V0)/alpha(1));
    z2 = @(x) 1/l0(2)*(x - a - (E(i) - V1)/alpha(2));
    
    % left end x = 0
    ML = [1, 1;
          1i*k, -1i*k];
      
    % match exponential-linear boundary at x = 0
    M1 = [airy(0,z1(0)), airy(2,z1(0));
         1/l0(1)*airy(1,z1(0)), 1/l0(1)*airy(3,z1(0))];
     
    % match linear-linear boundary at x = a
    M2 = [airy(0,z1(a)), airy(2,z1(a));
         1/l0(1)*airy(1,z1(a)), 1/l0(1)*airy(3,z1(a))];
     
    % match linear-linear boundary at x = a
    M3 = [airy(0,z2(a)), airy(2,z2(a));
         1/l0(2)*airy(1,z2(a)), 1/l0(2)*airy(3,z2(a))];
    
    % match linear-exponential boundary at x = b
    M4 = [airy(0,z2(b)), airy(2,z2(b));
         1/l0(2)*airy(1,z2(b)), 1/l0(2)*airy(3,z2(b))];
     
    % right end at x = b
    MR = [exp(1i*k*b), 0;
         1i*k*exp(1i*k*b), 0];
     
    %Mtot = inv(ML) * M1 * inv(M2) * M3 * inv(M3) * M4 * inv(M4) * MR;
    Mtot = (ML\M1)*(M2\M3)*(M3\M4)*(M4\MR);

    T(i) = 1/Mtot(1,1);
    T(i) = conj(T(i))*T(i);

end

subplot(211)
plot(E,T,'Linewidth',2)
xlabel('E / Hartree')
ylabel('Transmission')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')

subplot(212)
v1 = @(x) alpha(1)*x + V0;
v2 = @(x) alpha(2)*(x - a) + V1;
y1 = linspace(0,a,100); y2 = linspace(a,a+b,100);
plot([y1,y2],[v1(y1),v2(y2)])

%% Trapezoidal barrier v2 (1 line)

clear all
clc
close all

m = 100;

x_start = 100;

V0L = 10;

V0 = 1;
V1 = 10;
a = 1;

alpha = (V1 - V0)/a ;

E = linspace(1,50,1000) * V0;

T = zeros(1,length(E));
T_exact = zeros(1,length(E));
for i = 1:length(E)
    
    k = sqrt(2*m*(E(i)-V0L));
    k2 = conj(k)*k;
    
    l01 = (2*m*alpha)^(-1/3);
    z1 = @(x) 1/l01*(x - x_start - (E(i) - V0)/alpha);
    
    ML = [exp(1i*k*x_start), exp(-1i*k*x_start);
          1i*k*exp(1i*k*x_start), -1i*k*exp(-1i*k*x_start)];
      
    M1 = [airy(0,z1(x_start)), airy(2,z1(x_start));
         1/l01*airy(1,z1(x_start)), 1/l01*airy(3,z1(x_start))];
     
    M2 = [airy(0,z1(x_start+a)), airy(2,z1(x_start+a));
         1/l01*airy(1,z1(x_start+a)), 1/l01*airy(3,z1(x_start+a))];
     
    MR = [exp(1i*k*(x_start+a)), 0;
         1i*k*exp(1i*k*(x_start+a)), 0];
     
    %Mtot = inv(ML) * M1 * inv(M2) * MR;
    Mtot = (ML\M1)*(M2\MR);
   
    T(i) = 1/Mtot(1,1);
    T(i) = conj(T(i))*T(i);
    
    % exact result
    f = airy(0,z1(x_start));
    fp = 1/l01*airy(1,z1(x_start));
    g = airy(2,z1(x_start));
    gp = 1/l01*airy(3,z1(x_start));
    q = airy(0,z1(x_start+a));
    qp = 1/l01*airy(1,z1(x_start+a));
    r = airy(2,z1(x_start+a));
    rp = 1/l01*airy(3,z1(x_start+a));
    
    T_exact(i) = (2i*k*(r*qp - q*rp))/...
        ( (-k2*f*r - fp*rp + k2*g*q + gp*qp) + 1i*k*(fp*r - f*rp - gp*q + g*qp) );
    T_exact(i) = conj(T_exact(i))*T_exact(i);
    
end

yyaxis left
plot(E/V0,T_exact)
yyaxis right
plot(E/V0,T)

idx = find(~isnan(T));

fprintf('Error = %4.10f\n',norm(T(idx)-T_exact(idx)))


%% Triangular barrier (2 lines)

clear all
clc
close all

m = 100;
a = 1;
b = 1;

V0 = 1;
V1 = 1.1;

alpha = [(V1 - V0)/a, (V0 - V1)/b];

E = linspace(0,8,1000);

T = zeros(1,length(E));
for i = 1:length(E)
    
    k = sqrt(2*m*E(i));
    k2 = conj(k)*k;
    
    l0 = airylengthscale(alpha,m);
    z1 = @(x) 1/l0(1)*(x - (E(i) - V0)/alpha(1));
    z2 = @(x) 1/l0(2)*(x - a - (E(i) - V1)/alpha(2));
    
    % left end x = 0
    ML = [1, 1;
          1i*k, -1i*k];
      
    % match exponential-linear boundary at x = 0
    M1 = [airy(0,z1(0)), airy(2,z1(0));
         1/l0(1)*airy(1,z1(0)), 1/l0(1)*airy(3,z1(0))];
     
    % match linear-linear boundary at x = a
    M2 = [airy(0,z1(a)), airy(2,z1(a));
         1/l0(1)*airy(1,z1(a)), 1/l0(1)*airy(3,z1(a))];
     
    % match linear-linear boundary at x = a
    M3 = [airy(0,z2(a)), airy(2,z2(a));
         1/l0(2)*airy(1,z2(a)), 1/l0(2)*airy(3,z2(a))];
    
    % match linear-exponential boundary at x = b
    M4 = [airy(0,z2(b)), airy(2,z2(b));
         1/l0(2)*airy(1,z2(b)), 1/l0(2)*airy(3,z2(b))];
     
    % right end at x = b
    MR = [exp(1i*k*b), 0;
         1i*k*exp(1i*k*b), 0];
     
    %Mtot = inv(ML) * M1 * inv(M2) * M3 * inv(M3) * M4 * inv(M4) * MR;
    Mtot = (ML\M1)*(M2\M3)*(M3\M4)*(M4\MR);

    T(i) = 1/Mtot(1,1);
    T(i) = conj(T(i))*T(i);

end

subplot(211)
plot(E,T,'Linewidth',2)
xlabel('E / Hartree')
ylabel('Transmission')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')

subplot(212)
v1 = @(x) alpha(1)*x + V0;
v2 = @(x) alpha(2)*(x - a) + V1;
y1 = linspace(0,a,100); y2 = linspace(a,a+b,100);
plot([y1,y2],[v1(y1),v2(y2)])

%% Triangular barrier v2 (2 lines)

clear all
clc
close all

m = 100;
a = 0.5;
b = 1;


V0L = 0.5;
V0R = 1;

V0 = -2;
V1 = 1.5;

alpha = [(V1 - V0)/a, (V0 - V1)/b];

[M1] = @(x,e) airyMatrix(alpha(1),x,V0,0,e,m);
[M2] = @(x,e) airyMatrix(alpha(2),x,V1,a,e,m);

E = linspace(0,8,1000);

T = zeros(1,length(E));
for i = 1:length(E)
    
    e = E(i);
    
    kL = sqrt(2*m*(e-V0L));
    kR = sqrt(2*m*(e-V0R));

    % left end x = 0
    ML = [1, 1;
          1i*kL, -1i*kL];
     
    % right end at x = b
    MR = [exp(1i*kR*b), 0;
         1i*kR*exp(1i*kR*b), 0];
     
    % WRONG!
    %Mtot = (ML\M1(0,e))*(M1(a,e)\M2(a,e))*(M2(a,e)\M2(b,e))*(M2(b,e)\MR);
    
    % RIGHT?
    Mtot = (ML\M1(0,e))*(M1(a,e)\M2(a,e))*(M2(b,e)\MR);

    T(i) = 1/Mtot(1,1) * kR / kL;
    T(i) = conj(T(i))*T(i);

end

subplot(211)
plot(E,T,'Linewidth',2)
xlabel('E / Hartree')
ylabel('Transmission')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')

subplot(212)
v1 = @(x) alpha(1)*x + V0;
v2 = @(x) alpha(2)*(x - a) + V1;
y1 = linspace(0,a,100); y2 = linspace(a,a+b,100);
plot([y1,y2],[v1(y1),v2(y2)])

%% Triangular barrier v3 (2 lines)

clear all
clc
close all

m = 100;
a = 0.5;
b = 1;

x0 = 0;

V0L = 0.5;
V0R = 1;

V0 = 0.2;
V1 = 5;

alpha = [(V1 - V0)/a, (V0 - V1)/b];

[M1] = @(x,e) airyMatrix(alpha(1),x,V0,x0,e,m);
[M2] = @(x,e) airyMatrix(alpha(2),x,V1,x0+a,e,m);

E = linspace(0,8,1000);

T = zeros(1,length(E));
for i = 1:length(E)
    
    e = E(i);
    
    kL = sqrt(2*m*(e-V0L));
    kR = sqrt(2*m*(e-V0R));

    % left end x = 0
    ML = [exp(1i*kL*x0), exp(-1i*kL*x0);
          1i*kL*exp(1i*kL*x0), -1i*kL*exp(-1i*kL*x0)];
     
    % right end at x = b
    MR = [exp(1i*kR*(x0+b)), 0;
         1i*kR*exp(1i*kR*(x0+b)), 0];
     
    % WRONG!
    %Mtot = (ML\M1(0,e))*(M1(a,e)\M2(a,e))*(M2(a,e)\M2(b,e))*(M2(b,e)\MR);
    
    % RIGHT?
    Mtot = (ML\M1(x0,e))*(M1(x0+a,e)\M2(x0+a,e))*(M2(x0+b,e)\MR);

    T(i) = 1/Mtot(1,1) * kR / kL;
    T(i) = conj(T(i))*T(i);

end

subplot(211)
plot(E,T,'Linewidth',2)
xlabel('E / Hartree')
ylabel('Transmission')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')

subplot(212)
v1 = @(x) alpha(1)*x + V0;
v2 = @(x) alpha(2)*(x - a) + V1;
y1 = linspace(0,a,100); y2 = linspace(a,a+b,100);
plot([y1,y2],[v1(y1),v2(y2)])


%% General case
clear all
clc
close all

N = 40; 

x = linspace(-20,200,500);

m = 1;

V0L = 1;
V0R = 0;

%xrange = [0, 100];
%f = @(x) (3.15 - (x./100).^2) .* window(x,xrange(1),xrange(2));

% set = 0 at x = 0
% 0 = -a(0-p)^2 + b => -ap^2 + b = 0 => b = ap^2
a = 0.01;
x0 = 10;
b = a*x0^2;
xrange = [0, sqrt(b/a)+x0];
f = @(x) (-a*(x-x0).^2 + b) .* window(x,xrange(1),xrange(2));

% a = 0.01;
% b = 1;
% xrange = [-sqrt(b/a), sqrt(b/a)];
% f = @(x) (-a*x.^2 + b) .* window(x,xrange(1),xrange(2));

[x_interp, V_interp, Vdata] = PiecewiseLinearInterpolate(f,xrange,N);

x_interior = x_interp(2:end-1);

subplot(212)
hold on
for i = 1:length(Vdata)
    y = linspace(x_interp(i),x_interp(i+1),20);
    v1 = Vdata{i}.slope*y + Vdata{i}.intercept;
    plot(y,v1,'r-','Linewidth',2)
end
plot(x,f(x),'k-','Linewidth',1)
plot(x_interp,V_interp,'ro')
hold off
xlabel('x / a.u.')
ylabel('Potential / Hartree')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')
grid on


M = cell(1,N);
for i = 1:N  
    M{i} = @(x,e) airyMatrix(Vdata{i}.slope,x,Vdata{i}.ystart,Vdata{i}.xstart,e,m);
end

E = linspace(0,5,200);

T = zeros(1,length(E));
for i = 1:length(E)
    
    e = E(i);
    
    kL = sqrt(2*m*(e-V0L));
    kR = sqrt(2*m*(e-V0R));
    

    % left end x = x_interp(1)
    ML = [  exp(1i*kL*x_interp(1)),       exp(-1i*kL*x_interp(1));
          1i*kL*exp(1i*kL*x_interp(1)), -1i*kL*exp(-1i*kL*x_interp(1))];
     
    % right end at x = x_interp(end)
    MR = [exp(1i*kR*x_interp(end)-x_interp(1)),     0;
         1i*kR*exp(1i*kR*x_interp(end)-x_interp(1)), 0];
     
    Mtr = eye(2);
    for j = 1:N-1
        Mtr = Mtr * (M{j}(x_interior(j)-x_interp(1),e)\M{j+1}(x_interior(j)-x_interp(1),e));
    end
    Mtot = (ML\M{1}(x_interp(1),e)) * Mtr * (M{end}(x_interp(end)-x_interp(1),e)\MR);
    %Mtot = (ML\M1(0,e))*(M1(a,e)\M2(a,e))*(M2(b,e)\MR);

%     Mtr = (ML\M{1}(x_interp(1),e)) * ...
%           (M{1}(x_interp(2),e)\M{2}(x_interp(2),e)) * ...
%           (M{2}(x_interp(2),e)\M{2}(x_interp(3),e)) * ...
%           (M{2}(x_interp
    %Mtot = (ML\M1(0,e))*(M1(a,e)\M2(a,e))*(M2(a,e)\M2(b,e))*(M2(b,e)\MR);

    T(i) = 1/Mtot(1,1) * kR / kL;
    T(i) = conj(T(i))*T(i);

end

subplot(211)
semilogy(E,T,'Linewidth',2)
xlabel('E / Hartree')
ylabel('log_{10}(|T|^2) Transmission')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')
grid on


%%

function [M] = airyMatrix(alpha,x,y_start,x_start,e,m)

    l0 = airylengthscale(alpha,m);
    
    z = @(x) 1/l0*(x - x_start - (e - y_start)/alpha);
    
    M = [airy(0,z(x)), airy(2,z(x));
         1/l0*airy(1,z(x)), 1/l0*airy(3,z(x))];
     
end

function [l0] = airylengthscale(alpha,m)
    
    l0 = zeros(1,length(alpha));
    for i = 1:length(alpha)
        l0(i) = (2*m*abs(alpha(i))).^(-1/3);
        if alpha(i) < 0
            l0(i) = -l0(i);
        end
    end
    
end
    
function [W] = window(x,low,high)

    W = heaviside(x - low) .* heaviside(high - x);
    
end
    