function [ xmin, fmin, results ] = golden_section_v2( f, a, b )

% 1D minimization using the Golden Section Method
% Input:
%   f = (anonymous) function of one variable to be minimzed, must be
%        unimodal within the domain [a,b]
%   a = lower-end domain boundary
%   b = upper-end domain boundary


% Minimum interval length for stopping criterion
tol = 1e-4;             

% Initialize parameters
L = b-a;  
x1 = a + (b-a)*0.382;   
x2 = a + (b-a)*0.618;
f1 = f(x1);
f2 = f(x2);
niter = 0;
results = [];
results(1,:) = [0, a, b,f1,f2];

% Begin loop
while (L > tol)

    if f1 > f2        
    % Eliminate the region [a, x1] and update parameters
        a = x1;
        x1 = x2;
        f1 = f2;
        x2 = a + (b-a)*0.618;
        f2 = f(x2);
    else 
    % Eliminate the region [x2, b] and update parameters
        b = x2;
        x2 = x1;
        f2 = f1;
        x1 = a + (b-a)*0.382;
        f1 = f(x1);
    end
    
    % Calculate final interval length
    L = abs(b-a);
    
    % Update iteration counter
    niter = niter + 1;
    
    % Store the results
    results(niter+1,:) = [niter + 1, a, b, f1, f2];
    

end

% When small enough, calculate xmin as average of end interval
xmin = (a+b)/2;
fmin = f(xmin);

% Write optimization results to text file
% fileID = fopen('opthist_gs.txt','w');
% fprintf(fileID,'%6.0s %6.0s %6.0s %6.0s %6.0s\n','#', 'a', 'b', 'f1', 'f2');
% for i = 1:length(results(:,1))
%     fprintf(fileID,'%6.4f %6.4f %6.4f %6.4f %6.4f\n',results(i,:));
% end
% fclose(fileID);



end

