function [xout] = grad_vAB(A,B,atom_coordinates,atom_valency,RX)
    
    xout = zeros(3,1);
    
    for i = 1:length(A.exps)
        for j = 1:length(B.exps)
%             N1 = norm_factor(A.shell,A.exps(i));
%             N2 = norm_factor(B.shell,B.exps(j));
            temp = 0;
            for k = 1:size(atom_coordinates,1)
                temp = temp - atom_valency(k)*derivative_nuclear_attraction(A.exps(i),A.shell,A.origin,B.exps(j),B.shell,B.origin,...
                                                                            atom_coordinates(k,:),RX);
            end
            xout = xout + A.coef(i)*B.coef(j)*temp;
        end
    end
    

end