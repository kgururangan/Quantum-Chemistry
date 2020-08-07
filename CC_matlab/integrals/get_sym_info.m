function [G,mult_table] = get_sym_info(point_group)

    % identity matrix
    I = eye(3);

    % rotation matrix about z axis
    Rz = @(x) [cos(x),  sin(x), 0;
               -sin(x), cos(x), 0;
                 0,       0,    1];

    % generates matrix that reflects across plane with UNIT normal n
    mirror_plane = @(n) [1-2*n(1)^2, -2*n(1)*n(2), -2*n(1)*n(3);
                        -2*n(1)*n(2), 1-2*n(2)^2, -2*n(2)*n(3);
                        -2*n(3)*n(1), -2*n(2)*n(3), 1-2*n(3)^2];    


switch point_group

    case {'C2V','c2v'}
        
        G{1} = I;

        G{2} = Rz(pi);

        G{3} = mirror_plane([0,1,0]);

        G{4} = mirror_plane([1,0,0]);
        
        mult_table = [1, 2, 3, 4;
                      2, 1, 4, 3;
                      3, 4, 1, 2;
                      4, 3, 2, 1];
    otherwise
        disp('Enter valid point group!')
end

end

