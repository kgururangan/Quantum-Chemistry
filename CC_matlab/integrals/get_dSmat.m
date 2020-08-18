
function [dSmat_cell] = get_dSmat(basis_cell,grad_coords)

    Nat = size(grad_coords,1);
    
    dSmat_cell = cell(Nat,1);
    
    for M = 1:Nat

        dSmat = zeros(length(basis_cell),length(basis_cell),3);
        
        RX = grad_coords(M,:);

        for I = 1:length(basis_cell)
            for J = I:length(basis_cell)
                dSmat(I,J,:) = grad_sAB(basis_cell{I},basis_cell{J},RX);
                dSmat(J,I,:) = dSmat(I,J,:);
            end
        end
        
        dSmat_cell{M} = dSmat;
        
    end
end

