function [] = query_orbital_space(sys_cc)

    Norb = sys_cc.Norb/2; % total number spatial orbitals
    Ncore = sys_cc.Ncore/2; % number of core spatial orbitals
    Nfzc = sys_cc.NFZ_core/2; % number of frozen occupied spatial orbitals
    Nfzv = sys_cc.NFZ_vir/2; % number of frozen unoccupied spatial orbitals
    Nact_h = sys_cc.Nact_h/2; % number of active occupied spatial orbitals
    Nact_p = sys_cc.Nact_p/2; % number of active unoccupied spatial orbitals
    Nelec = sys_cc.Nelec;
    
    
    fprintf('Frozen Core: ')
    for i = 1:Nfzc
        if i < Nfzc
            fprintf(' %d, ',i);
        else
            fprintf(' %d ',i);
        end
    end
    fprintf('\n')
    fprintf('Core: ')
    for i  = Nfzc+1:Ncore+Nfzc
        if i < Ncore+Nfzc
            fprintf(' %d, ',i);
        else
            fprintf(' %d ',i);
        end
    end
    for i = Ncore+Nfzc+1:Nelec

end

