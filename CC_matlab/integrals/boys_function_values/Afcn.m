function [A] = Afcn(n,alpha)

    if n >= alpha
        
        A = 0.0;
        
        kmax = 5;
        
        for k = 0:kmax
            A
            A = A + alpha^k/factorial(n+k+1);
        end
        
        A = A*factorial(n)*exp(-alpha);
        
        

        
    else
        
        A = (factorial(n) - igamma(n+1,alpha))/alpha^(n+1);
        
    end
    
end