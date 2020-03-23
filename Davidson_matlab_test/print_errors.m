function [ ] = print_errors( V_true, W_true, e_true, V, W, e )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    fprintf('\nPrinting errors in eigenpairs...\n\n')
    for j = 1:length(e_true)

        ov_right = (V(:,j)'*V_true(:,j));
        theta_right = acosd(ov_right);

        ov_left = (W(:,j)'*W_true(:,j));
        theta_left = acosd(ov_left);


        err_th_right = min([abs(theta_right-180),abs(theta_right)]);
        err_th_left = min([abs(theta_left-180),abs(theta_left)]);

        err_e = abs(e_true(j) - e(j));

        fprintf('Root-%d\n',j)
        fprintf('Error in eigenvalue = %4.10f   Right Vec rotation = %4.2f deg   Left Vec rotation = %4.2f deg\n',err_e, theta_right, theta_left)

    end


end

