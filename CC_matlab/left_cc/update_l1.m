function [L1] = update_l1(L1, L2, HBar, FM, VM, omega, occ, unocc, flag_ground)

    Nocc = length(occ); Nunocc = length(unocc);
    
    [X_ia, ~] = build_L_HBar(L1, L2, HBar, FM, VM, flag_ground, 1);
    
    
    for i = 1:Nocc
        for a = 1:Nunocc
            MP = HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i);
            L1(i,a) = -X_ia(i,a)/(MP - omega);    
        end
    end
    
end

