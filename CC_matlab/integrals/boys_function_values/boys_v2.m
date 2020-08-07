function [val] = boys_v2(m,z)


    N = 5;
    Q = 36;
    
    u = 1.0;
    a = u/N;
    
    val = 0.0;
    
    for s = 0:N-1
        
        temp = 0.0;
        
        for k = 0:Q
            
            Aval = Afcn(2*m+2*k,2*z*s*a^2);
            
            temp = temp + (power(-z,2*k)*power(a,2*m+2*k+1)* Aval)/factorial(k);
                    
        end
        
        val = val + temp*exp(-s^2*a^2);
        
    end



end

function [A] = Afcn(n,alpha)

    if n > alpha
        
        A = 0.0;
        
        kmax = 30;
        
        for k = 0:kmax
            A = A + alpha^k/factorial(n+k+1);
        end
        
        A = A*factorial(n)*exp(-alpha);
        
    else
        
        A = (factorial(n) - igamma(n+1,alpha))/alpha^(n+1);
        
    end
    
end





