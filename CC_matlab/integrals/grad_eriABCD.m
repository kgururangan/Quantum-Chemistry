function [xout] = grad_eriABCD(A,B,C,D,RX)
    
    xout = zeros(3,1);
    
    for i = 1:length(A.exps)
        for j = 1:length(B.exps)
            for k = 1:length(C.exps)
                for l = 1:length(D.exps)
%                     N1 = norm_factor(A.shell,A.exps(i));
%                     N2 = norm_factor(B.shell,B.exps(j));
%                     N3 = norm_factor(C.shell,C.exps(k));
%                     N4 = norm_factor(D.shell,D.exps(l));
                    
                    xout = xout + ...
                        A.coef(i)*B.coef(j)*C.coef(k)*D.coef(l)*...
                        derivative_electron_repulsion(A.exps(i),A.shell,A.origin,B.exps(j),B.shell,B.origin,C.exps(k),C.shell,C.origin,D.exps(l),D.shell,D.origin,RX);
                end
            end
        end
    end
    
end

