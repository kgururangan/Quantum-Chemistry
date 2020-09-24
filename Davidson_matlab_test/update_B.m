function [add_B] = update_B(r,B,add_B,thresh_vec)

    r_orth = ortho_root_vec(r,B(:,1:curr_size));
    if norm(r_orth)/norm(r) > thresh_vec
        B(:,curr_size+1) = r_orth/norm(r_orth);
        curr_size = curr_size + 1;
    end


end

