function [orbs, Rat, idx_gnd, idx_exc] = build_orbitals(num_chromophore, xi_gnd, cm_gnd, xi_exc, cm_exc, origin, chromophore_spacing, flag_exc)

    Rat = zeros(num_chromophore,3);
    
    if flag_exc == 1
        orbs = cell(1,2*num_chromophore);
        ct = 1;
        for i = 1:num_chromophore
            Rat(i,:) = origin+(i-1)*chromophore_spacing*[1,0,0];
            g.origin = origin+(i-1)*chromophore_spacing*[1,0,0];
            g.exps = xi_gnd{i};
            g.coeff = cm_gnd{i};
            g.coeff = g.coeff*normalize_AO(g);
            orbs{ct} = g;
            ct = ct + 1;
        end

        for i = 1:num_chromophore
            g.origin = origin+(i-1)*chromophore_spacing*[1,0,0];
            g.exps = xi_exc{i};
            g.coeff = cm_exc{i};
            g.coeff = g.coeff*normalize_AO(g);
            orbs{ct} = g;
            ct = ct + 1;
        end

        idx_gnd = 1:num_chromophore;
        idx_exc = num_chromophore+1:2*num_chromophore;
    else
        orbs = cell(1,num_chromophore);
        ct = 1;
        for i = 1:num_chromophore
            Rat(i,:) = origin+(i-1)*chromophore_spacing*[1,0,0];
            g.origin = origin+(i-1)*chromophore_spacing*[1,0,0];
            g.exps = xi_gnd{i};
            g.coeff = cm_gnd{i};
            g.coeff = g.coeff*normalize_AO(g);
            orbs{ct} = g;
            ct = ct + 1;
        end
        idx_gnd = 1:num_chromophore;
        idx_exc = [];
    end
end



