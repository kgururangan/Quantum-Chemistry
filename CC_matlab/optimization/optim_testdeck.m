clear all
clc
close all

a = 20;
b = 100;

% rosenbrock function
f = @(x) (a-x(1))^2 + b*(x(2)-x(1)^2)^2;
df = @(x) [-2*(a-x(1)) - 4*b*x(1)*(x(2)-x(1)^2);
           2*b*(x(2)-x(1)^2)];
       
      

x0 = rand(2,1);

opts.tol = 1e-6;
opts.maxit = 15e3;
opts.alpha = [0, 2];
H0 = eye(size(x0,1));

[xsol] = bfgs(f,df,x0,H0,opts);
%[xmin, fmin] = conjugate_gradient(f, x0)

%%
clear all
clc
close all

ff = @(x) fcn1(x)

a = ff([1,2])

function [cell_out] = fcn1(x)
    cell_out{1} = x(1).^2;
    cell_out{2} = 5*x(2);
end