% function [val_large,val_small] = boys_v3(n,x)
% 
% 
%     kmax = 6;
%     
%     val_large = fact2(2*n-1)/power(2,n+1)*sqrt(pi/power(x,2*n+1));
%     val_small = 0.0;
%     for k = 0:kmax
%         val_small = power(-x,k)/(factorial(k)*(2*n+2*k+1));
%     end
% 
% end

function [val] = boys_analytical(n,x)

    if x == 0 
        val = 1/(2*n+1);
        return
    end
        
    PHI = @(t) normcdf(t,0,1);

    t1 = sqrt(pi)*power(2,-n)*prod([1:2:2*n-1])*(2*PHI(sqrt(2*x))-1);
    t2 = 0;
    for k = 0:n-1
        m = n-k-1;
        if m == 0
            t2 = t2 + power(-x,k);
        end
        if m > 0
            t2 = t2 + prod([0.5-n,-k-1.5])*power(-x,k);
        end
    end
    val = 0.5*power(x,-n-0.5)*(t1+t2*power(-1,n)*sqrt(x)*exp(-x));

end

% t1=sqrt(pi)*2**(-n)*prod(seq(1,2*n-1,2))*(2*pnorm(sqrt(2*x))-1) 
% t2=0
% for (k in seq(0,n-1,1))
% {m=n-k-1
% if (m==0) t2=t2+(—x)**k
% if (m>0) t2=t2+prod(seq(0.5-n,-k-1.5))*(-x)**k} 
% tt=0.5*x**(-n-0.5)*(t1+t2*(-1)**n*sqrt(x)*exp(—x)) return(tt)
