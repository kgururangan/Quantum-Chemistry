function [c] = sym_product(MO_idx1, MO_idx2, irreps, sym_group)

    if strcmp(sym_group,'c1')
        group_mult_table = ones(2);
    elseif strcmp(sym_group,'c2v')
        group_mult_table = [1, 2, 3, 4;
                            2, 1, 4, 3
                            3, 4, 1, 2
                            4, 3, 2, 1];
    else
        disp('Symmetry %s not supported',sym_group)
    end
    
    
    c = group_mult_table(irreps(MO_idx1),irreps(MO_idx2));

end

