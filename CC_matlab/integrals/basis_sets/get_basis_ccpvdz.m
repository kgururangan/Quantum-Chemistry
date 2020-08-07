function [basis_cell] = get_basis_ccpvdz(element_cell, atom_coordinates)

    basis_cell = {};

    for J = 1:length(element_cell)
        
        element_label = element_cell{J};


        switch element_label

            case {'H','HYDROGEN','Hydrogen'}

                atom_cell = cell(1,5);

                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [1.301000E+01, 1.962000E+00, 4.446000E-01, 1.220000E-01];
                atom_cell{1}.coef = [1.968500E-02, 1.379770E-01, 4.781480E-01, 5.012400E-01];

                atom_cell{2}.shell = [0, 0, 0];
                atom_cell{2}.exps = [1.220000E-01];
                atom_cell{2}.coef = [1.000000E+00];

                atom_cell{3}.shell = [1, 0, 0];
                atom_cell{3}.exps = [7.270000E-01];
                atom_cell{3}.coef = [1.000000E+00];
                atom_cell{4}.shell = [0, 1, 0];
                atom_cell{4}.exps = [7.270000E-01];
                atom_cell{4}.coef = [1.000000E+00];
                atom_cell{5}.shell = [0, 0, 1];
                atom_cell{5}.exps = [7.270000E-01];
                atom_cell{5}.coef = [1.000000E+00];

            case {'O','OXYGEN','Oxygen'}

                atom_cell = cell(1,15);

                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [1.172000E+04, 1.759000E+03, 4.008000E+02, 1.137000E+02, 3.703000E+01, 1.327000E+01, ...
                                      5.025000E+00, 1.013000E+00, 3.023000E-01];
                atom_cell{1}.coef = [7.100000E-04, 5.470000E-03, 2.783700E-02, 1.048000E-01, 2.830620E-01, 4.487190E-01, ...
                                      2.709520E-01, 1.545800E-02, -2.585000E-03];

                atom_cell{2}.shell = [0, 0, 0];
                atom_cell{2}.exps = [1.172000E+04, 1.759000E+03, 4.008000E+02, 1.137000E+02, 3.703000E+01, 1.327000E+01, ...
                                      5.025000E+00, 1.013000E+00, 3.023000E-01];
                atom_cell{2}.coef = [-1.600000E-04, -1.263000E-03, -6.267000E-03, -2.571600E-02, -7.092400E-02, -1.654110E-01, ...
                                      -1.169550E-01, 5.573680E-01, 5.727590E-01];

                atom_cell{3}.shell = [0, 0, 0];
                atom_cell{3}.exps = [3.023000E-01];
                atom_cell{3}.coef = [1.000000E+00];

                atom_cell{4}.shell = [1, 0, 0];
                atom_cell{4}.exps = [1.770000E+01, 3.854000E+00, 1.046000E+00, 2.753000E-01];
                atom_cell{4}.coef = [4.301800E-02, 2.289130E-01, 5.087280E-01, 4.605310E-01];
                atom_cell{5}.shell = [0, 1, 0];
                atom_cell{5}.exps = [1.770000E+01, 3.854000E+00, 1.046000E+00, 2.753000E-01];
                atom_cell{5}.coef = [4.301800E-02, 2.289130E-01, 5.087280E-01, 4.605310E-01];
                atom_cell{6}.shell = [0, 0, 1];
                atom_cell{6}.exps = [1.770000E+01, 3.854000E+00, 1.046000E+00, 2.753000E-01];
                atom_cell{6}.coef = [4.301800E-02, 2.289130E-01, 5.087280E-01, 4.605310E-01];

                atom_cell{7}.shell = [1, 0, 0];
                atom_cell{7}.exps = [2.753000E-01];
                atom_cell{7}.coef = [1.000000E+00];
                atom_cell{8}.shell = [0, 1, 0];
                atom_cell{8}.exps = [2.753000E-01];
                atom_cell{8}.coef = [1.000000E+00];
                atom_cell{9}.shell = [0, 0, 1];
                atom_cell{9}.exps = [2.753000E-01];
                atom_cell{9}.coef = [1.000000E+00];

                atom_cell{10}.shell = [2, 0, 0];
                atom_cell{10}.exps = [1.185000E+00];
                atom_cell{10}.coef = [1.0000000];
                atom_cell{11}.shell = [0, 2, 0];
                atom_cell{11}.exps = [1.185000E+00];
                atom_cell{11}.coef = [1.0000000];
                atom_cell{12}.shell = [0, 0, 2];
                atom_cell{12}.exps = [1.185000E+00];
                atom_cell{12}.coef = [1.0000000];
                atom_cell{13}.shell = [1, 1, 0];
                atom_cell{13}.exps = [1.185000E+00];
                atom_cell{13}.coef = [1.0000000];
                atom_cell{14}.shell = [0, 1, 1];
                atom_cell{14}.exps = [1.185000E+00];
                atom_cell{14}.coef = [1.0000000];
                atom_cell{15}.shell = [1, 0, 1];
                atom_cell{15}.exps = [1.185000E+00];
                atom_cell{15}.coef = [1.0000000];
                
            case{'C','CARBON','Carbon'}
                atom_cell = cell(1,15);

                atom_cell{1}.shell = [0, 0, 0];
                atom_cell{1}.exps = [6.665000E+03, 1.000000E+03, 2.280000E+02, 6.471000E+01, 2.106000E+01, 7.495000E+00, ...
                                      2.797000E+00, 5.215000E-01, 1.596000E-01];
                atom_cell{1}.coef = [6.920000E-04, 5.329000E-03, 2.707700E-02, 1.017180E-01, 2.747400E-01, 4.485640E-01, ...
                                      2.850740E-01, 1.520400E-02, -3.191000E-03];

                atom_cell{2}.shell = [0, 0, 0];
                atom_cell{2}.exps = [6.665000E+03, 1.000000E+03, 2.280000E+02, 6.471000E+01, 2.106000E+01, 7.495000E+00, ...
                                      2.797000E+00, 5.215000E-01, 1.596000E-01];
                atom_cell{2}.coef = [-1.460000E-04, -1.154000E-03, -5.725000E-03, -2.331200E-02, -6.395500E-02, -1.499810E-01, ...
                                      -1.272620E-01, 5.445290E-01, 5.804960E-01];

                atom_cell{3}.shell = [0, 0, 0];
                atom_cell{3}.exps = [1.596000E-01];
                atom_cell{3}.coef = [1.000000E+00];

                atom_cell{4}.shell = [1, 0, 0];
                atom_cell{4}.exps = [9.439000E+00, 2.002000E+00 , 5.456000E-01 , 1.517000E-01];
                atom_cell{4}.coef = [3.810900E-02, 2.094800E-01, 5.085570E-01, 4.688420E-01];
                atom_cell{5}.shell = [0, 1, 0];
                atom_cell{5}.exps = [9.439000E+00, 2.002000E+00 , 5.456000E-01 , 1.517000E-01];
                atom_cell{5}.coef = [3.810900E-02, 2.094800E-01, 5.085570E-01, 4.688420E-01];
                atom_cell{6}.shell = [0, 0, 1];
                atom_cell{6}.exps = [9.439000E+00, 2.002000E+00 , 5.456000E-01 , 1.517000E-01];
                atom_cell{6}.coef = [3.810900E-02, 2.094800E-01, 5.085570E-01, 4.688420E-01];

                atom_cell{7}.shell = [1, 0, 0];
                atom_cell{7}.exps = [1.517000E-01];
                atom_cell{7}.coef = [1.000000E+00];
                atom_cell{8}.shell = [0, 1, 0];
                atom_cell{8}.exps = [1.517000E-01];
                atom_cell{8}.coef = [1.000000E+00];
                atom_cell{9}.shell = [0, 0, 1];
                atom_cell{9}.exps = [1.517000E-01];
                atom_cell{9}.coef = [1.000000E+00];

                atom_cell{10}.shell = [2, 0, 0];
                atom_cell{10}.exps = [5.500000E-01];
                atom_cell{10}.coef = [1.0000000];
                atom_cell{11}.shell = [0, 2, 0];
                atom_cell{11}.exps = [5.500000E-01];
                atom_cell{11}.coef = [1.0000000];
                atom_cell{12}.shell = [0, 0, 2];
                atom_cell{12}.exps = [5.500000E-01];
                atom_cell{12}.coef = [1.0000000];
                atom_cell{13}.shell = [1, 1, 0];
                atom_cell{13}.exps = [5.500000E-01];
                atom_cell{13}.coef = [1.0000000];
                atom_cell{14}.shell = [0, 1, 1];
                atom_cell{14}.exps = [5.500000E-01];
                atom_cell{14}.coef = [1.0000000];
                atom_cell{15}.shell = [1, 0, 1];
                atom_cell{15}.exps = [5.500000E-01];
                atom_cell{15}.coef = [1.0000000];

            otherwise 

                fprintf('\nElement %s does not have cc-pVDZ supported!',EL)
        end
        
        for I = 1:size(atom_cell,2)
            atom_cell{I}.origin = atom_coordinates(J,:);
        end
        
        basis_cell = cat(2,basis_cell,atom_cell);
        
    end
            
end

