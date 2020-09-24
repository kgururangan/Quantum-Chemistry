function [basis_cell] = get_basis_sto3G(element_cell, atom_coordinates)

    basis_cell = {};

    for J = 1:length(element_cell)
        
        element_label = element_cell{J};


        switch element_label

            case {'H','HYDROGEN','Hydrogen'}

                atom_cell = cell(1,1);

                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00];
                atom_cell{1}.coef = [0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00];


            case {'O','OXYGEN','Oxygen'}

                atom_cell = cell(1,5);

                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [0.1307093214E+03, 0.2380886605E+02, 0.6443608313E+01];
                atom_cell{1}.coef = [0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00];

                atom_cell{2}.shell = [0, 0, 0];
                atom_cell{2}.exps = [0.5033151319E+01,  0.1169596125E+01, 0.3803889600E+00];
                atom_cell{2}.coef = [-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00];
                
                atom_cell{3}.shell = [1, 0, 0];
                atom_cell{3}.exps = [0.5033151319E+01,  0.1169596125E+01, 0.3803889600E+00];
                atom_cell{3}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];
                atom_cell{4}.shell = [0, 1, 0];
                atom_cell{4}.exps = [0.5033151319E+01,  0.1169596125E+01, 0.3803889600E+00];
                atom_cell{4}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];
                atom_cell{5}.shell = [0, 0, 1];
                atom_cell{5}.exps = [0.5033151319E+01,  0.1169596125E+01, 0.3803889600E+00];
                atom_cell{5}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];


                
            case{'C','CARBON','Carbon'}
                
                atom_cell = cell(1,5);
                
                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [0.7161683735E+02, 0.1304509632E+02, 0.3530512160E+01];
                atom_cell{1}.coef = [0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00];
                                 
                atom_cell{2}.shell = [0, 0, 0];
                atom_cell{2}.exps = [0.2941249355E+01, 0.6834830964E+00, 0.2222899159E+00];
                atom_cell{2}.coef = [-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00];
                
                atom_cell{3}.shell = [1, 0, 0];
                atom_cell{3}.exps = [0.2941249355E+01, 0.6834830964E+00, 0.2222899159E+00];
                atom_cell{3}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];
                atom_cell{4}.shell = [0, 1, 0];
                atom_cell{4}.exps = [0.2941249355E+01, 0.6834830964E+00, 0.2222899159E+00];
                atom_cell{4}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];
                atom_cell{5}.shell = [0, 0, 1];
                atom_cell{5}.exps = [0.2941249355E+01, 0.6834830964E+00, 0.2222899159E+00];
                atom_cell{5}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];
       
            case{'Si','SILICON','Silicon'}
                
                atom_cell = cell(1,9);
                
                atom_cell{1}.shell = [0,0,0];
                atom_cell{1}.exps = [0.4077975514E+03, 0.7428083305E+02, 0.2010329229E+02];
                atom_cell{1}.coef = [0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00];
                
                atom_cell{2}.shell = [0,0,0];
                atom_cell{2}.exps = [0.2319365606E+02, 0.5389706871E+01, 0.1752899952E+01];
                atom_cell{2}.coef = [-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00];
                atom_cell{3}.shell = [1,0,0];
                atom_cell{3}.exps = [0.2319365606E+02, 0.5389706871E+01, 0.1752899952E+01];
                atom_cell{3}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];
                atom_cell{4}.shell = [0,1,0];
                atom_cell{4}.exps = [0.2319365606E+02, 0.5389706871E+01, 0.1752899952E+01];
                atom_cell{4}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];
                atom_cell{5}.shell = [0,0,1];
                atom_cell{5}.exps = [0.2319365606E+02, 0.5389706871E+01, 0.1752899952E+01];
                atom_cell{5}.coef = [0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00];

                atom_cell{6}.shell = [0,0,0];
                atom_cell{6}.exps = [0.1478740622E+01, 0.4125648801E+00, 0.1614750979E+00];
                atom_cell{6}.coef = [-0.2196203690E+00, 0.2255954336E+00, 0.9003984260E+00];
                atom_cell{7}.shell = [1,0,0];
                atom_cell{7}.exps = [0.1478740622E+01, 0.4125648801E+00, 0.1614750979E+00];
                atom_cell{7}.coef = [0.1058760429E-01, 0.5951670053E+00, 0.4620010120E+00];
                atom_cell{8}.shell = [0,1,0];
                atom_cell{8}.exps = [0.1478740622E+01, 0.4125648801E+00, 0.1614750979E+00];
                atom_cell{8}.coef = [0.1058760429E-01, 0.5951670053E+00, 0.4620010120E+00];
                atom_cell{9}.shell = [0,0,1];
                atom_cell{9}.exps = [0.1478740622E+01, 0.4125648801E+00, 0.1614750979E+00];
                atom_cell{9}.coef = [0.1058760429E-01, 0.5951670053E+00, 0.4620010120E+00];

            otherwise 

                fprintf('\nElement %s does not have sto-3G supported!',EL)
        end
        
        for I = 1:size(atom_cell,2)
            atom_cell{I}.origin = atom_coordinates(J,:);
        end
        
        basis_cell = cat(2,basis_cell,atom_cell);
        
    end
            
end