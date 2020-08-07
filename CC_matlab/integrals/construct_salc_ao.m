function [] = construct_salc_ao(basis,point_group)

    [sym_group, mult_table] = get_sym_info(point_group);
    H = length(sym_group);
    Norb = length(basis);

    CN = zeros(H*Norb,Norb);

    % loop over each AO in the basis
    for i = 1:Norb
        orb = basis{i};
        for j = 1:H
            G = sym_group{j};
            
        


end

function [sign] = get_orbital_sym_phase(orb,G)
    % return the phase factor characterizing the various
    % ao's like s, p, and d. We do so using the fact that the
    % angular variant portion of these function is in the 
    % cartesian prefactor. i.e.
    % s -> invariant to all rotation (sign = +1 always)
    % px -> x  py -> y  pz -> z
    % dx2 -> x^2 dy2 -> y^2  dz2 -> z^2  dxy -> xy  dyz -> yz  dxz -> xz


    if orb.shell == [0,0,0]
        sign = 1;
    end




end