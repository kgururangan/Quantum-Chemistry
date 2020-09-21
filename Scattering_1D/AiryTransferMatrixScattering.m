clear all
clc
close all

N = 40; 

x = linspace(-20,200,500);

m = 1;

V0L = 0.1; % can only have scattering states for E > V0L
V0R = 0.1; % can only have scattering states for E > V0R

% xrange = [0, 100];
% f = @(x) (3.15 - (x./100).^2) .* window(x,xrange(1),xrange(2));

% set = 0 at x = 0
% 0 = -a(0-p)^2 + b => -ap^2 + b = 0 => b = ap^2
a = 0.002;
x0 = 10;
b = a*x0^2;
xrange = [0, sqrt(b/a)+x0];
f = @(x) (-a*(x-x0).^2 + b) .* window(x,xrange(1),xrange(2));

% a = 0.001;
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

E = linspace(max(V0R,V0L),1,500);

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
%         TMP = (M{j}(x_interior(j)-x_interp(1),e)\M{j+1}(x_interior(j)-x_interp(1),e));
%         err = norm(M{j}(x_interior(j)-x_interp(1),e)*TMP - M{j+1}(x_interior(j)-x_interp(1),e));
%         if err > 1e-10
%             break;
%         end
       % Mtr = Mtr * (pinv(M{j}(x_interior(j)-x_interp(1),e))*M{j+1}(x_interior(j)-x_interp(1),e));
    end
    Mtot = (ML\M{1}(x_interp(1),e)) * Mtr * (M{end}(x_interp(end)-x_interp(1),e)\MR);

    T(i) = 1/Mtot(1,1) * kR / kL;
    T(i) = conj(T(i))*T(i);

end

subplot(211)
%semilogy(E,T,'Linewidth',2)
plot(E,T,'Linewidth',2)
xlabel('E / Hartree')
ylabel('log_{10}(|T|^2) Transmission')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')
grid on

%% Functions

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
    