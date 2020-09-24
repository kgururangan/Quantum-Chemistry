function [xmin, fmin] = conjugate_gradient(f, x0, grad_tol)

% Unconstrained multivariable minimization using Fletcher-Reeves conjugate
% gradient method

% Input: f = @(x,y,z) ... anonymous function of 2 or 3 variables
%        x0 = [a,b,c]' ... 3x1 vector of initial guess position
%        grad_tol = gradient stopping condition (default 1e-4)
% Output: xmin ... 3x1 vector of minimum position
%         fmin ... value of function at minimum
%         results ... niter x 10 matrix, 
%         [Iteration #, x(1),x(2),x(3),f(x),d(1),d(2),d(3), alpha, |grad(f)|]

fprintf('\n=====================++Entering CG Routine++=====================\n')



% Initialize parameters
x = x0;
h = 1e-4;          % Derivative step size
niter = 0;
results = [];
nvar = size(x0,1);

if nargin < 3
    grad_tol = 1e-4;
end

% nvar x nvar identity matrix
Id = eye(nvar);

% Store initial value of gradient and its magnitude
ck_old = zeros(nvar,1); ck = zeros(nvar,1);
for j = 1:nvar
    ck_old(j) = 1/h*(f(x+h/2*Id(:,j))-f(x-h/2*Id(:,j)));
end    
ck_mag = sqrt(dot(ck_old,ck_old));

% Store initial value of search direction
d_old = -ck_old;

%results(1,:) = [0, x0(1), x0(2), x0(3), f(x0(1),x0(2),x0(3)), d_old(1), d_old(2), d_old(3), NaN, ck_mag];

while ck_mag > grad_tol
    
    % Calculate gradient using central difference formula
    for j = 1:nvar
        ck(j) = 1/h*(f(x+h/2*Id(:,j))-f(x-h/2*Id(:,j)));
    end    
    
    ck_mag = sqrt(dot(ck,ck));
    ck_old_mag = sqrt(dot(ck_old,ck_old));
    
    % Fletcher-Reeves beta factor for updating search direction
    beta = (ck_mag/ck_old_mag)^2;
    
    % Update search direction
    d = -ck + beta*d_old;
    
    % Optimize one-variable function of step size
    fa = @(a) f(x + a*d);
    [alpha, fmin, ~] = golden_section_v2(fa, 0, 10);
    
    % Update minimum position
    x = x + alpha*d;
    
    % Store value of ck and d for next iteration
    ck_old = ck;
    d_old = d;
    
    % Update iteration counter
    niter = niter + 1;
    
    % Store the results
    %results(niter+1,:) = [niter,x(1),x(2),x(3),fmin,d(1),d(2),d(3),alpha,ck_mag];
    fprintf('Iter - %d:      FVAL = %4.8f      OPTIM = %4.8f\n',niter, fmin, ck_mag)
    
end

xmin = x;

if ck_mag <= grad_tol
    fprintf('CG successfuly converged in %d iterations\n',niter-1)
else
    fprintf('CG not converged\n')
end

% Write optimization results to text file
% fileID = fopen('opthist.txt','w');
% fprintf(fileID,'%6.2s %6.2s %6.2s %6.2s %6.2s %6.2s %6.2s %6.2s %6.5s %6.2s\n','#','x1','x2','x3','f','d1','d2','d3','alpha','df');
% for i = 1:length(results(:,1))
%     fprintf(fileID,'%6.0f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',results(i,:));
% end
% fclose(fileID);


end
