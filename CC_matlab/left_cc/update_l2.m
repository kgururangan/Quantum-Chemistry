function [L2] = update_l2(L1, L2, HBar, FM, VM, omega, occ, unocc, flag_ground)

    Nocc = length(occ); Nunocc = length(unocc);
    
    [~, X_ijab] = build_L_HBar(L1, L2, HBar, FM, VM, flag_ground, 1);
    
    
    for i = 1:Nocc
        for j = i+1:Nocc
            for a = 1:Nunocc
                for b = a+1:Nunocc
                    
                    MP =   HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i)...
                          +HBar{1}{2,2}(b,b) - HBar{1}{1,1}(j,j);
                
                    L2(i,j,a,b) = -X_ijab(i,j,a,b)./(MP - omega);
                    L2(j,i,a,b) = -L2(i,j,a,b);
                    L2(i,j,b,a) = -L2(i,j,a,b);
                    L2(j,i,b,a) = L2(i,j,a,b);
                    
                end
            end
        end
    end
    
end

