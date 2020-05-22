function [MO_sp] = spatial_to_spinorb(MO)

    Norb = size(MO,1);
    
    if length(size(MO)) == 2
        MO_sp = zeros(2*Norb,2*Norb);
        for i = 1:2*Norb
            for j = 1:2*Norb
                if mod(i,2) == mod(j,2)  
                    if mod(i,2) == 1
                        i0 = floor(i/2)+1;
                    else
                        i0 = floor(i/2);
                    end
                    if mod(j,2) == 1
                        j0 = floor(j/2)+1;
                    else
                        j0 = floor(j/2);
                    end
                    MO_sp(i,j) = MO(i0,j0);
                end
            end
        end
    
            
    else
        MO_sp = zeros(2*Norb,2*Norb,2*Norb,2*Norb);
        for i = 1:2*Norb
            for j = 1:2*Norb
                for k = 1:2*Norb
                    for l = 1:2*Norb
                        if mod(i,2) == mod(k,2) && mod(j,2) == mod(l,2)
                            if mod(i,2) == 1
                                i0 = floor(i/2)+1;
                            else
                                i0 = floor(i/2);
                            end
                            if mod(j,2) == 1
                                j0 = floor(j/2)+1;
                            else
                                j0 = floor(j/2);
                            end
                            if mod(k,2) == 1
                                k0 = floor(k/2)+1;
                            else
                                k0 = floor(k/2);
                            end
                            if mod(l,2) == 1
                                l0 = floor(l/2)+1;
                            else
                                l0 = floor(l/2);
                            end
                            MO_sp(i,j,k,l) = MO(i0,j0,k0,l0);
                        end
                    end
                end
            end
        end
    end
    
end
