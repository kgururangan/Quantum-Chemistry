function [dets_char] = parse_determinants(dets, num_chromophore)

    n_elec = 2*num_chromophore;
    norb = 2*n_elec;
    
    dets_char = zeros(size(dets,1),1); % 1 if S, 2 if CT, 3 if TT
    
    idx_alpha = [1:2:norb];
    idx_beta = [2:2:norb];
    idx_at = zeros(num_chromophore, 4);
    
    ct = 1;
    for i = 1:2:n_elec
        idx_at(ct,:) = [i, i+1, i+n_elec, i+n_elec+1];
        ct = ct+1;
    end
    
    for i = 1:size(dets,1)
        D = dets(i,:);
        temp = zeros(num_chromophore,3);
        for j = 1:num_chromophore
            [vals_at,~] = intersect(D,idx_at(j,:),'stable');
            [vals_alpha,~] = intersect(vals_at,idx_alpha,'stable');
            [vals_beta,~] = intersect(vals_at,idx_beta,'stable');
            [vals_exc_at,~] = intersect(vals_at,[n_elec+1:norb],'stable');
            temp(j,1) = length(vals_at);
            temp(j,2) = length(vals_alpha);
            temp(j,3) = length(vals_beta);
            temp(j,4) = length(vals_exc_at);
        end
        
        dets_char(i) = assign_character(temp);
        

    end

end

function [character] = assign_character(mat)
    % added to remove high energy determinants
    if any(mat(:,1) == 0) || sum(mat(:,4)) >= 3 || sum(mat(:,4)) == 0
        character = 0;
    elseif any(mat(:,1) == 3)
        character = 2;
    elseif any(mat(:,2)==2) || any(mat(:,3)==2)
        character = 3;
    elseif mat(1,2:end) == mat(2,2:end)
        character = 1;
    else
        character = 0;
    end

end
