function [xout] = tAB(A,B)
    
    xout = 0.0;
    
    for i = 1:length(A.exps)
        for j = 1:length(B.exps)
            N1 = norm_factor(A.shell,A.exps(i));
            N2 = norm_factor(B.shell,B.exps(j));
            xout = xout + ...
                   N1*N2*A.coef(i)*B.coef(j)*kinetic(A.exps(i),A.shell,A.origin,B.exps(j),B.shell,B.origin);
        end
    end
    

end
