function [GmnE] = myFranckCondon(Bv, Dpv, n_modes, n_quanta)

S = 1;
for k = 1:n_modes
    
    B = Bv(k);
    dp = Dpv(k);
    
    G00E = sqrt(2*B/(1+B^2))*exp(-B^2*dp^2/(1+B^2)/2); 
    
    GmnEmt = zeros(n_quanta,n_quanta);
    
    x2 = -dp*B/sqrt((1-B^4));
    x2i = -1i*B*x2;
    
    for i = 0:n_quanta-1 % G states
        for j = 0:n_quanta-1 % E states
            a1 = ((1-B^2)/2/(1+B^2))^((i+j)/2);
            a2 = sqrt(factorial(j)*factorial(i));
            a3 = 1;
            X = 0;
            for l = 0:min(i,j)
                b1 = (4*B/(1-B^2))^l;
                b2 = (-1i)^(j-l)/factorial(l)/factorial(j-l)/factorial(i-l);
                X = X + b1*b2*myHermite(j-l,x2i)*myHermite(i-l,x2);
            end
            GmnEmt(i + 1,j + 1) = a1*a2*a3*X;
        end
    end
    GmnEm = GmnEmt*G00E; %matrix of F.C. factors
    S = kron(GmnEm,S);
end

GmnE = S;
EnmG = GmnE';

end