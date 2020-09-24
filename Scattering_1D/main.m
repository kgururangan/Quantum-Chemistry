clear all
clc
close all

N = 40;

m = 10;

V0L = 0; % can only have scattering states for E > V0L
V0R = 1; % can only have scattering states for E > V0R

% bound states
% xrange = [0, 100];
% f = @(x) 2.15*(1 - exp(-x./20));

xrange = [0, 100];
f = @(x) (3.15 - (x./100).^2);

% set = 0 at x = 0
% 0 = -a(0-p)^2 + b => -ap^2 + b = 0 => b = ap^2
% a = 0.001;
% x0 = 20;
% b = a*x0^2;
% xrange = [0, sqrt(b/a)+x0];
% f = @(x) (-a*(x-x0).^2 + b);% .* window(x,xrange(1),xrange(2));

% Eckart potential
% xrange = [0,50]; x = linspace(xrange(1),xrange(2),500);
% a = 100; s0 = 0;
% mat = [exp(2*pi*(xrange(1)-s0)/a)/(1+exp(2*pi*(xrange(1)-s0)/a)), exp(2*pi*(xrange(1)-s0)/a)/(1+exp(2*pi*(xrange(1)-s0)/a))^2;
%        exp(2*pi*(xrange(2)-s0)/a)/(1+exp(2*pi*(xrange(2)-s0)/a)), exp(2*pi*(xrange(2)-s0)/a)/(1+exp(2*pi*(xrange(2)-s0)/a))^2];
% b = [V0L+0.0, V0R-2]';
% Vpar = mat\b;
% f = @(x) Vpar(1)*exp(2*pi*x/a)./(1+exp(2*pi*x/a)) + Vpar(2)*exp(2*pi*x/a)./(1+exp(2*pi*x/a)).^2;

%
%
E = linspace(0,5,500);

[T] = AiryTMM(f,xrange,E,m,V0L,V0R,N);




