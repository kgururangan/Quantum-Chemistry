function [Zmat] = get_Zmat(basis_cell,atom_coordinates,atom_valency)

    Zmat = zeros(length(basis_cell));
    
    for I = 1:length(basis_cell)
        for J = I:length(basis_cell)
            Zmat(I,J) = tAB(basis_cell{I},basis_cell{J}) + vAB(basis_cell{I},basis_cell{J},atom_coordinates,atom_valency);
            Zmat(J,I) = Zmat(I,J);
        end
    end
    
end