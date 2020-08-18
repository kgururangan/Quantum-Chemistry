function [xout] = grad_sAB(A,B,RX)
    
    xout = zeros(3,1);
    
    for i = 1:length(A.exps)
        for j = 1:length(B.exps)
%             N1 = norm_factor(A.shell,A.exps(i));
%             N2 = norm_factor(B.shell,B.exps(j));
            xout = xout + ...
                   A.coef(i)*B.coef(j)*derivative_overlap(A.exps(i),A.shell,A.origin,B.exps(j),B.shell,B.origin,RX);
        end
    end
    

end

