function [hf_energy] = calculate_hf_energy(e1int,e2int,Nocc_a,Nocc_b)

    % Calculate the scf energy
    hf_energy = 0.0;
    for i = 1:Nocc_a
        hf_energy = hf_energy + e1int(i,i);
    end
    for i = 1:Nocc_b
        hf_energy = hf_energy + e1int(i,i);
    end
    for i = 1:Nocc_a
        for j = 1:Nocc_b
            hf_energy = hf_energy + e2int(i,j,i,j);
        end
    end
    for i = 1:Nocc_a
        for j = 1:Nocc_a
            hf_energy = hf_energy + 0.5*(e2int(i,j,i,j)-e2int(i,j,j,i));
        end
    end
    for i = 1:Nocc_b
        for j = 1:Nocc_b
            hf_energy = hf_energy + 0.5*(e2int(i,j,i,j)-e2int(i,j,j,i));
        end
    end

end

