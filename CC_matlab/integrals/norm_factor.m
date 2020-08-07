function [N] = norm_factor(shell,alpha)

    l = shell(1); m = shell(2); n = shell(3);
    
    N = sqrt(power(2,2*(l+m+n)+1.5)*power(alpha,l+m+n+1.5)...
             /fact2(2*l-1)/fact2(2*m-1)/fact2(2*n-1)/power(pi,1.5));
     
end