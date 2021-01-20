function [f] = boys(n,T)

    if n == 0 % deal with F0(T) separately since we have analytic form
        if T > 0
            f = 0.5*sqrt(pi/T)*erf(sqrt(T));
        else
            f = 1;
        end
  
    else % n > 0

        % numerical implementation - FAST :)
        f = hyp1F1(n+0.5,n+1.5,-T)/(2*n+1);
        %f = boys_analytical(n,T);


        if isnan(f) ||  f > 1/(2*n+1)
            
            % Analytical Boys function from 
            %fprintf('Using analytical Boys...\n')
            f = boys_analytical(n,T);
            %f = hyp1F1(n+0.5,n+1.5,-T)/(2*n+1);
 
        
            if isnan(f) || f > 1/(2*n+1)
                % native Matlab function - symbolic and SLOW!
                % pFq(a,b,z) = hypergeom(a, b, z) where length(a) = p and length(b) = q
                fprintf('Boys function NaN - using symbolic implementation...\n')
                f = hypergeom(n+0.5,n+1.5,-T)/(2*n+1);
            end
            
        end
        
    end


end


% def boys(n, x):
%     if x > 0.0:
%         f = 2.0*x**(n + 0.5)
%         g = gamma(n + 0.5)
%         gi = 1.0 - gammainc(n + 0.5, x, regularized=True)
%         return g*gi/f
%     else:
%         return 1.0/(n*2 + 1)