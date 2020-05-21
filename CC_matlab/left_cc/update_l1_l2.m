function [L1,L2] = update_l1_l2(X_ia,X_ijab,omega,HBar,occ,unocc)

% updates using abij order

    Nocc = length(occ); Nunocc = length(unocc);
    
    L1 = zeros(Nunocc,Nocc); L2 = zeros(Nunocc,Nunocc,Nocc,Nocc);
    
    for i = 1:Nocc
        for a = 1:Nunocc
            
            MP = HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i);

            L1(a,i) = -X_ia(a,i)./(MP - omega);

            for j = i+1:Nocc
                for b = a+1:Nunocc
                    
                    MP =   HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i)...
                          +HBar{1}{2,2}(b,b) - HBar{1}{1,1}(j,j);
                      
                    L2(a,b,i,j) = -X_ijab(a,b,i,j)./(MP - omega);                   
                    L2(b,a,i,j) = -L2(a,b,i,j);
                    L2(a,b,j,i) = -L2(a,b,i,j);
                    L2(b,a,j,i) = L2(a,b,i,j);
                end
            end
        end
    end

end
