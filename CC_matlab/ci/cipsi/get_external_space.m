function [V_ext] = get_external_space(Wavefunction,sys)

    V_ext = [];
    for K = 1:size(Wavefunction.Dets,1)
        D = Wavefunction.Dets(K,:);
        [singles,doubles] = get_singles_and_doubles(D,sys);
        V_ext = [V_ext; singles.dets; doubles.dets];
    end
    
    ct = 1;
    idx_repeat = [];
    for i = 1:size(V_ext,1)
        for j = 1:size(V_ext,2)
            if i ~= j
                if V_ext(i,:) == V_ext(j,:)
                    idx_repeat(ct) = j;
                    ct = ct + 1;                   
                end
            end
        end
    end
   
    if ~isempty(idx_repeat)
        V_ext(idx_repeat,:) = [];
    end
                
                
end

