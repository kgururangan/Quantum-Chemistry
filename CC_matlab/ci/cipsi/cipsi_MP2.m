function [eMP2_array] = cipsi_MP2(Wavefunction,V_ext,sys)

    eMP2_array = zeros(size(V_ext,1),1);
    
    for K = 1:size(V_ext,1)
        [temp] = calculate_mp2(V_ext(K,:),Wavefunction,sys);
        eMP2_array(K) = temp;
    end
    
end

