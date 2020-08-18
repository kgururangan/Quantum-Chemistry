function [basis_cell] = get_basis_631G(element_cell, atom_coordinates)

    basis_cell = {};

    for J = 1:length(element_cell)
        
        element_label = element_cell{J};


        switch element_label

            case {'H','HYDROGEN','Hydrogen'}

                atom_cell = cell(1,2);

                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [0.1873113696E+02, 0.2825394365E+01, 0.6401216923E+00];
                atom_cell{1}.coef = [0.3349460434E-01, 0.2347269535E+00, 0.8137573261E+00];

                atom_cell{2}.shell = [0, 0, 0];
                atom_cell{2}.exps = [0.1612777588E+00];
                atom_cell{2}.coef = [1.000000E+00];

            case {'O','OXYGEN','Oxygen'}

                atom_cell = cell(1,9);

                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [0.5484671660E+04, 0.8252349460E+03, 0.1880469580E+03, ...
                                     0.5296450000E+02, 0.1689757040E+02, 0.5799635340E+01];
                atom_cell{1}.coef = [0.1831074430E-02, 0.1395017220E-01, 0.6844507810E-01, ...
                                     0.2327143360E+00, 0.4701928980E+00, 0.3585208530E+00];

                atom_cell{2}.shell = [0, 0, 0];
                atom_cell{2}.exps = [0.1553961625E+02, 0.3599933586E+01, 0.1013761750E+01];
                atom_cell{2}.coef = [-0.1107775495E+00, -0.1480262627E+00, 0.1130767015E+01];
                
                atom_cell{3}.shell = [1, 0, 0];
                atom_cell{3}.exps = [0.1553961625E+02, 0.3599933586E+01, 0.1013761750E+01];
                atom_cell{3}.coef = [0.7087426823E-01, 0.3397528391E+00, 0.7271585773E+00];
                atom_cell{4}.shell = [0, 1, 0];
                atom_cell{4}.exps = [0.1553961625E+02, 0.3599933586E+01, 0.1013761750E+01];
                atom_cell{4}.coef = [0.7087426823E-01, 0.3397528391E+00, 0.7271585773E+00];
                atom_cell{5}.shell = [0, 0, 1];
                atom_cell{5}.exps = [0.1553961625E+02, 0.3599933586E+01, 0.1013761750E+01];
                atom_cell{5}.coef = [0.7087426823E-01, 0.3397528391E+00, 0.7271585773E+00];

                atom_cell{6}.shell = [0, 0, 0];
                atom_cell{6}.exps = [0.2700058226E+00];
                atom_cell{6}.coef = [1.000000E+00];

                atom_cell{7}.shell = [1, 0, 0];
                atom_cell{7}.exps = [0.2700058226E+00];
                atom_cell{7}.coef = [1.000000E+00];
                atom_cell{8}.shell = [0, 1, 0];
                atom_cell{8}.exps = [0.2700058226E+00];
                atom_cell{8}.coef = [1.000000E+00];
                atom_cell{9}.shell = [0, 0, 1];
                atom_cell{9}.exps = [0.2700058226E+00];
                atom_cell{9}.coef = [1.000000E+00];

                
            case{'C','CARBON','Carbon'}
                
                atom_cell = cell(1,9);
                
                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [0.3047524880E+04, 0.4573695180E+03, 0.1039486850E+03, ...
                                     0.2921015530E+02, 0.9286662960E+01, 0.3163926960E+01];
                atom_cell{1}.coef = [0.1834737132E-02, 0.1403732281E-01, 0.6884262226E-01, ...
                                     0.2321844432E+00, 0.4679413484E+00, 0.3623119853E+00];
                                 
                atom_cell{2}.shell = [0, 0, 0];
                atom_cell{2}.exps = [0.7868272350E+01, 0.1881288540E+01, 0.5442492580E+00];
                atom_cell{2}.coef = [-0.1193324198E+00, -0.1608541517E+00, 0.1143456438E+01];
                
                atom_cell{3}.shell = [1, 0, 0];
                atom_cell{3}.exps = [0.7868272350E+01, 0.1881288540E+01, 0.5442492580E+00];
                atom_cell{3}.coef = [0.6899906659E-01, 0.3164239610E+00, 0.7443082909E+00];
                atom_cell{4}.shell = [0, 1, 0];
                atom_cell{4}.exps = [0.7868272350E+01, 0.1881288540E+01, 0.5442492580E+00];
                atom_cell{4}.coef = [0.6899906659E-01, 0.3164239610E+00, 0.7443082909E+00];
                atom_cell{5}.shell = [0, 0, 1];
                atom_cell{5}.exps = [0.7868272350E+01, 0.1881288540E+01, 0.5442492580E+00];
                atom_cell{5}.coef = [0.6899906659E-01, 0.3164239610E+00, 0.7443082909E+00];
        
                atom_cell{6}.shell = [0, 0, 0];
                atom_cell{6}.exps = [0.1687144782E+00];
                atom_cell{6}.coef = [0.1000000000E+01];
                
                atom_cell{7}.shell = [1, 0, 0];
                atom_cell{7}.exps = [0.1687144782E+00];
                atom_cell{7}.coef = [0.1000000000E+01];
                atom_cell{8}.shell = [0, 1, 0];
                atom_cell{8}.exps = [0.1687144782E+00];
                atom_cell{8}.coef = [0.1000000000E+01];
                atom_cell{9}.shell = [0, 0, 1];
                atom_cell{9}.exps = [0.1687144782E+00];
                atom_cell{9}.coef = [0.1000000000E+01];

            otherwise 

                fprintf('\nElement %s does not have 6-31G supported!',EL)
        end
        
        for I = 1:size(atom_cell,2)
            atom_cell{I}.origin = atom_coordinates(J,:);
        end
        
        basis_cell = cat(2,basis_cell,atom_cell);
        
    end
            
end
