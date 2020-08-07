function [Smat] = get_Smat(basis_cell)

    Smat = zeros(length(basis_cell));
    
    for I = 1:length(basis_cell)
        for J = I:length(basis_cell)
            Smat(I,J) = sAB(basis_cell{I},basis_cell{J});
            Smat(J,I) = Smat(I,J);
        end
    end
end

