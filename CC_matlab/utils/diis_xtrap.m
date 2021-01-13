function [X_xtrap] = diis_xtrap(X_list, diis_resid_list)

        [vec_dim, diis_dim] = size(X_list);
        
        B_dim = diis_dim + 1;
        
        B = zeros(B_dim);
        B(end,:) = -1;
        B(:,end) = -1;
        B(end,end) = 0;
        
        for i = 1:diis_dim
            for j = i:diis_dim
                % residual inner product using flattened vectors
                B(i,j) = diis_resid_list(:,i)'*diis_resid_list(:,j);
                B(j,i) = B(i,j);
            end
        end
        
        rhs = zeros(B_dim,1);
        rhs(end) = -1;
        
        %if rcond(B) > eps % THIS RUINS THINGS! SKIPPING DIIS IS BAD!
            %coeff = B\rhs;
            coeff = solveGauss(B,rhs);
            X_xtrap = zeros(vec_dim,1);
            for k = 1:diis_dim
                X_xtrap = X_xtrap + coeff(k)*X_list(:,k);
            end
        %else
        %    X_xtrap = X_list(:,end);
        %end
        
end
            
