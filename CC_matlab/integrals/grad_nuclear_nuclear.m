function [dVnn_cell] = grad_nuclear_nuclear(atom_coordinates,atom_valency,grad_coords)

    Nat = size(atom_coordinates,1);
    Nat_grad = size(grad_coords,1);
    dVnn_cell = cell(Nat_grad,1);
    
    for M = 1:Nat_grad
        
        dVnn = zeros(3,1);
        
        for A = 1:Nat
            if all(atom_coordinates(A,:) == grad_coords(M,:))
                for B = 1:Nat
                    if A ~= B
                        Ra = atom_coordinates(A,:)';
                        Rb = atom_coordinates(B,:)';
                        Rab = norm(Ra-Rb);
                        dVnn = dVnn - atom_valency(A)*atom_valency(B)*...
                            (Ra-Rb)/Rab^3;
                    end
                end
            end
        end
        
        dVnn_cell{M} = dVnn;
        
    end
    
end
