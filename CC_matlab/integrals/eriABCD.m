function [xout] = eriABCD(A,B,C,D)
    
    xout = 0.0;
    
    for i = 1:length(A.exps)
        for j = 1:length(B.exps)
            for k = 1:length(C.exps)
                for l = 1:length(D.exps)
                    N1 = norm_factor(A.shell,A.exps(i));
                    N2 = norm_factor(B.shell,B.exps(j));
                    N3 = norm_factor(C.shell,C.exps(k));
                    N4 = norm_factor(D.shell,D.exps(l));
                    
                    xout = xout + ...
                        N1*N2*N3*N4*...
                        A.coef(i)*B.coef(j)*C.coef(k)*D.coef(l)*...
                        electron_repulsion(A.exps(i),A.shell,A.origin,B.exps(j),B.shell,B.origin,C.exps(k),C.shell,C.origin,D.exps(l),D.shell,D.origin);
                end
            end
        end
    end
    
end

