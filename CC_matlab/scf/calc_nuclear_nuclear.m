function [Vnuc] = calc_nuclear_nuclear(atom_coordinates,Z)

%     atom_coordinates = molecule.atom_coordinates;
%     Z = molecule.atom_valency;

    Vnuc = 0.0;
    Nat = length(Z);
    
    for i = 1:Nat
        for j = 1:Nat
            if i ~= j
                Ri = atom_coordinates(i,:);
                Rj = atom_coordinates(j,:);
                Vnuc = Vnuc + Z(i)*Z(j)/norm(Ri-Rj);
            end
        end
    end

    Vnuc = 0.5*Vnuc;
    
end
    

