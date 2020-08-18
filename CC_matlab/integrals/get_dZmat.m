function [dZmat_cell] = get_dZmat(basis_cell,atom_coordinates,atom_valency,grad_coords)

    Nat = size(grad_coords,1);
    
    dZmat_cell = cell(Nat,1);
    
    for M = 1:Nat

        dZmat = zeros(length(basis_cell),length(basis_cell),3);
        
        RX = grad_coords(M,:);

        for I = 1:length(basis_cell)
            for J = I:length(basis_cell)
                dZmat(I,J,:) = grad_tAB(basis_cell{I},basis_cell{J},RX) + grad_vAB(basis_cell{I},basis_cell{J},atom_coordinates,atom_valency,RX);
                dZmat(J,I,:) = dZmat(I,J,:);
            end
        end
        
        dZmat_cell{M} = dZmat;
        
    end
    
end