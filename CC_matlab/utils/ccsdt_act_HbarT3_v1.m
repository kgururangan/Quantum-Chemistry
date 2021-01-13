% function [] = ccsdt_act_HBarT3_A_writer(out_proj)

%     % (HBar*T3)_C    
%     D1 = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
%     D2 = einsum_kg(H1A.vv,t3a,'ce,abeijk->abcijk');
%     D3 = 0.5*einsum_kg(H2A.oooo,t3a,'mnij,abcmnk->abcijk'); 
%     D4 = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
%     D5 = einsum_kg(H2A.voov,t3a,'cmke,abeijm->abcijk');
%     D6 = einsum_kg(H2B.voov,t3b,'cmke,abeijm->abcijk');
%     
%     % A(k/ij)
%     D13 = D1 + D3;
%     D13 = D13 - permute(D13,[1,2,3,6,5,4]) - permute(D13,[1,2,3,4,6,5]);
%     
%     % A(c/ab)
%     D24 = D2 + D4;
%     D24 = D24 - permute(D24,[3,2,1,4,5,6]) - permute(D24,[1,3,2,4,5,6]);
%    
%     % A(k/ij)A(c/ab)
%     D56 = D5 + D6;
%     D56 = D56 - permute(D56,[1,2,3,6,5,4]) - permute(D56,[1,2,3,4,6,5])...
%               - permute(D56,[3,2,1,4,5,6]) - permute(D56,[1,3,2,4,5,6]) ...
%               + permute(D56,[3,2,1,6,5,4]) + permute(D56,[3,2,1,4,6,5]) ...
%               + permute(D56,[1,3,2,6,5,4]) + permute(D56,[1,3,2,4,6,5]);    
%      
%      X3A_abcijk = MM23A + D13 + D24 + D56;
    
            
        %  D1 = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
        %  D1 = D1 - permute(D1,[1,2,3,6,5,4]) - permute(D1,[1,2,3,4,6,5]);

        for a = ['a','A']
            for b = ['b','B']
                for c = ['c','C']
                    for i = ['i','I']
                        for j = ['j','J']
                            for k = ['k','K']
                                
                                fprintf('\n%++++++++++%s%s%s%s%s%s++++++++++\n',a,b,c,i,j,k)
                                out_proj = {a,b,c,i,j,k};
                                d1 = write_HBarT3_oo(out_proj);
                                write_term = @(x,y) write_einsum(x,cell2mat(out_proj),y);
                                write_term(d1,'D1')
                                
                            end
                        end
                    end
                end
            end
        end
                                        
            

%end




function [d1] = write_HBarT3_A(out_proj)

    d1 = {};
    
    a = out_proj{1}; b = out_proj{2}; c = out_proj{3}; 
    i = out_proj{4}; j = out_proj{5}; k = out_proj{6};
    
    % H1(oo) * T3(vvvooo)
    % spin-integrated:
    %  D1 = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
    %  D1 = D1 - permute(D1,[1,2,3,6,5,4]) - permute(D1,[1,2,3,4,6,5]);
    %  m = act/unact
    
    ct = 1;
    arr1 = [m,k];
    arr2 = [a,b,c,i,j,m];
    
    ik = get_act_idx(k); ii = get_act_idx(i); ij = get_act_idx(j);
    if ik == ii && ik == ij % A(k/ij)
        arr1 = [m,k; m,i; m,j];
    elseif ik == ii && ik ~= ij
        arr1 = [m,k; m,i];
    elseif ik == ij && ik ~= ii
        arr1 = [m,k; m,j];
    else
        arr1 = [m,k];
    end
    
    for np = 1:size(arr1,1)
        idx1 = arr1(np,:);

        for m = ['m','M']
            coef = '';
            q1 = 'h1a';
            q2 = 't3a';
            arr1 = [m,k];
            arr2 = [a,b,c,i,j,m];
            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
            term1 = [sign,coef,q1,'(',arr1,')'];
            term2 = [q2,'(',new_arr,')'];
            d1{ct} = [term1,',',term2];
            ct = ct + 1;
        end
    end


end

function [new_arr,sign_char] = fix_t3_indices(arr,spin)

    sign = 1;

    switch spin
        
        case {'A','a'}
            part = arr(1:3); hole = arr(4:6);
            [perm_part] = fix_particle_idx(part);
            [perm_hole] = fix_hole_idx(hole);
            new_arr = arr([perm_part,perm_hole]);
            sign = sign * permutationParity(perm_part) * permutationParity(perm_hole);
            
        case {'B','b'}
            part = arr(1:2); hole = arr(4:5);
            [perm_part] = fix_particle_idx(part);
            [perm_hole] = fix_hole_idx(hole);
            new_arr = arr([perm_part,3,perm_hole,6]);
            sign = sign * permutationParity(perm_part) * permutationParity(perm_hole); 
            
        case {'c','C'}
            part = arr(2:3); hole = arr(5:6);
            [perm_part] = fix_particle_idx(part);
            [perm_hole] = fix_hole_idx(hole);
            new_arr = arr([1,perm_part+1,4,perm_hole+1]); % !!!
            sign = sign * permutationParity(perm_part) * permutationParity(perm_hole);

        case {'D','d'}
            part = arr(1:3); hole = arr(4:6);
            [perm_part] = fix_particle_idx(part);
            [perm_hole] = fix_hole_idx(hole);
            new_arr = arr([perm_part,perm_hole]);
            sign = sign * permutationParity(perm_part) * permutationParity(perm_hole);
            
        otherwise 
            disp('Spin type not allowed for T3!')
            
    end
    
    if sign == 1
        sign_char = '';
    else
        sign_char = '-';
    end
            

end

function [perm] = fix_particle_idx(idx)

    act_idx = zeros(1,length(idx));
    for i = 1:length(idx)
        val = get_act_idx(idx(i));
        act_idx(i) = val;
    end
    
    [~,perm] = sort(act_idx,'descend');

end

function [perm] = fix_hole_idx(idx)

    act_idx = zeros(1,length(idx));
    for i = 1:length(idx)
        val = get_act_idx(idx(i));
        act_idx(i) = val;
    end
    
    [~,perm] = sort(act_idx,'ascend');
    perm = perm + 3; %!!!
 
end

function [val] = get_character(char)

    val_act = get_act_idx(char);
    val_ph = get_ph_idx(char);
    
    if contains(char,'~') && val_ph == 1 % all occupied
        val = 5;
    elseif contains(char,'~') && val_ph == 0 % all unoccupied
        val = 6;
    else
    
        if val_act == 1 && val_ph == 1 % active particle
            val = 1;
        end

        if val_act == 1 && val_ph == 0 % active hole
            val = 2;
        end

        if val_act == 0 && val_ph == 1 % virtual
            val = 3;
        end

        if val_act == 0 && val_ph == 0 % core
            val = 4;
        end
    end
        
end

function [val] = get_act_idx(char)

    if length(char) > 1
        char1 = char(1);
    else
        char1 = char;
    end
    
    if isstrprop(char1,'upper')
        val = 1;
    else
        val = 0;
    end
       
end

function [val] = get_ph_idx(char)


    if length(char) > 1
        char1 = char(1);
    else
        char1 = char;
    end

    switch char1
        
        case {'m','n','i','j','k','M','N','I','J','K'}
            val = 0;
            
        case {'e','f','a','b','c','E','F','A','B','C'}
            val = 1;
            
    end
              
end

function [name] = get_subname(idx_cell)

    name = '';
    slices = cell(1,length(idx_cell));
    for i = 1:length(idx_cell)
        
        val = get_character(idx_cell{i});
        %ph = get_ph_idx(idx_cell{i});
        switch val
            case 1
                add_char = 'P';
                %name = [name,'P'];
            case 2
                add_char = 'H';
                %name = [name,'H'];
            case 3
                add_char = 'p';
                %name = [name,'p'];
            case 4
                add_char = 'h';
                %name = [name,'h'];
            case 5
                add_char = 'v';
                %name = [name,'v'];
            case 6
                add_char = 'o';
                %name = [name,'o'];
            otherwise
                fprintf('index character %s not recognized!\n',idx_cell{i})
        end
        
        name = [name,add_char];
        

    end

end

function [arr1, arr2, sign] = antisymm_arrs(arr_in_1, arr_in_2)
    
    n1 = length(arr_in_1);
    n2 = length(arr_in_2);
    nonsense = 100;
    temp1 = zeros(1,n1);
    temp2 = zeros(1,n2);
    for j = 1:n1
        val = get_character(arr_in_1(j));
        temp1(j) = val;
    end
    for j = 1:n2
        val = get_character(arr_in_2(j));
        temp2(j) = val;
    end
    if n1 < n2
        temp1 = [temp1, nonsense*ones(1,n2-n1)];
    end
    [stemp1,ix1] = sort(temp1,'ascend');
    [stemp2,ix2] = sort(temp2,'ascend');
    imatch = find(stemp1==stemp2);
    nmatch = length(imatch);
    im1 = imatch(ix1);
    im2 = imatch(ix2);
    arr1 = zeros(nmatch+1,n1); arr1(1,:) = arr_in_1;
    arr2 = zeros(nmatch+1,n2); arr2(1,:) = arr_in_2;
    sign = ones(1,nmatch);
    % im1 must be exchanged with im2 with proper sign of exchange
    for j = 1:nmatch
        temp1 = arr_in_1; 
        temp2 = arr_in_2;
        temp1(im1(j)) = temp2(im2(j));
        temp2(im2(j)) = arr_in_1(im1(j));      
        arr1(j+1,:) = temp1;
        arr2(j+1,:) = temp2;
    end
end
    
        
    
        





