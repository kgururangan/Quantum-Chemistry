function [singles,doubles] = get_singles_and_doubles(Det,sys)

    % note: sys.singles_dim and sys.doubles_dim enumerates ALL elements in
    % tensors C1 and C2. The determinants are ordered so you only need
    % permutationally unique ones. These will be generally less than (or
    % equal to) sys.singles_dim and sys.doubles_dim

    
    spinorbs = [1:2*sys.Norb];
    
    occ_alpha = Det(1:2:end-1); unocc_alpha = setdiff(spinorbs(1:2:end-1),occ_alpha);
    occ_beta = Det(2:2:end); unocc_beta = setdiff(spinorbs(2:2:end),occ_beta);
    
    % Singles
    ct_singles = 1;
    singles_dets = zeros(sys.singles_dim,sys.Nelec);
    singles_exc = zeros(sys.singles_dim,2);
    
    % C1A
    for i = 1:sys.Nocc_a
        for a = 1:sys.Nvir_a
            temp = Det;
            temp(occ_alpha(i)) = unocc_alpha(a);
            singles_dets(ct_singles,:) = sort(temp,'ascend');
            singles_exc(ct_singles,:) = [i,a];
            ct_singles = ct_singles + 1;
        end
    end
    
    % C1B
    for j = 1:sys.Nocc_b
        for b = 1:sys.Nvir_b
            temp = Det;
            temp(occ_beta(j)) = unocc_beta(b);
            singles_dets(ct_singles,:) = sort(temp,'ascend');
            singles_exc(ct_singles,:) = [j,b];
            ct_singles = ct_singles + 1;
        end
    end
    
    % Doubles
    ct_doubles = 1;
    doubles_dets = zeros(sys.doubles_dim,sys.Nelec);
    doubles_exc = zeros(sys.doubles_dim,4);
    
    %C2A
    for i = 1:sys.Nocc_a
        for j = i+1:sys.Nocc_a
            for a = 1:sys.Nvir_a
                for b = a+1:sys.Nvir_a
                    temp = Det;
                    temp(occ_alpha(i)) = unocc_alpha(a);
                    temp(occ_alpha(j)) = unocc_alpha(j);
                    doubles_dets(ct_doubles,:) = sort(temp,'ascend');
                    doubles_exc(ct_doubles,:) = [i, j, a, b];
                    ct_doubles = ct_doubles + 1;
                end
            end
        end
    end
    
    %C2B
    for i = 1:sys.Nocc_a
        for j = 1:sys.Nocc_b
            for a = 1:sys.Nvir_a
                for b = 1:sys.Nvir_b
                    temp = Det;
                    temp(occ_alpha(i)) = unocc_alpha(a);
                    temp(occ_beta(j)) = unocc_beta(j);
                    doubles_dets(ct_doubles,:) = sort(temp,'ascend');
                    doubles_exc(ct_doubles,:) = [i, j, a, b];
                    ct_doubles = ct_doubles + 1;
                end
            end
        end
    end
    
    %C2C
    for i = 1:sys.Nocc_b
        for j = i+1:sys.Nocc_b
            for a = 1:sys.Nvir_b
                for b = a+1:sys.Nvir_b
                    temp = Det;
                    temp(occ_beta(i)) = unocc_beta(a);
                    temp(occ_beta(j)) = unocc_beta(j);
                    doubles_dets(ct_doubles,:) = sort(temp,'ascend');
                    doubles_exc(ct_doubles,:) = [i, j, a, b];
                    ct_doubles = ct_doubles + 1;
                end
            end
        end
    end

            
    singles.dets = singles_dets(1:ct_singles-1,:); singles.exc = singles_exc(1:ct_singles-1,:);
    doubles.dets = doubles_dets(1:ct_doubles-1,:); doubles.exc = doubles_exc(1:ct_doubles-1,:);
    

%     occ = Det;
%     unocc = setdiff([1:sys.Norb],occ);
%     
%     % singles
%     ct_singles = 1;
%     for i = 1:length(occ)
%         for a = 1:length(unocc)
%             temp = Det;
%             if mod(occ(i),2) == mod(unocc(a),2)
%                 temp(find(temp==occ(i))) = unocc(a);
%                 exc_list_singles(ct_singles,:) = sort(temp,'ascend');
%                 ct_singles = ct_singles + 1;
%             end
%         end
%     end
%     exc_list_singles = unique(exc_list_singles,'rows');
%     
%     % doubles
%     ct_doubles = 1;
%     combs_occ = nchoosek(occ,2);
%     combs_unocc = nchoosek(unocc,2);
% 
%     
%     for I = 1:size(combs_occ,1)
%         i = combs_occ(I,1); j = combs_occ(I,2);
%         for A = 1:size(combs_unocc,1)
%             a = combs_unocc(A,1); b = combs_unocc(A,2);
%             if mod(i,2) == mod(a,2) && mod(j,2) == mod(b,2)
%                 temp = Det;
%                 temp(find(temp==occ(i))) = a; 
%                 temp(find(temp==occ(j))) = b;
%                 exc_list_doubles(ct_doubles,:) = sort(temp,'ascend');
%                 ct_doubles = ct_doubles + 1;
%             end
%             if mod(i,2) == mod(b,2) && mod(j,2) == mod(a,2)
%                 temp = Det;
%                 temp(find(temp==occ(i))) = b; 
%                 temp(find(temp==occ(j))) = a;
%                 exc_list_doubles(ct_doubles,:) = sort(temp,'ascend');
%                 ct_doubles = ct_doubles + 1;
%             end
%         end
%     end
%     exc_list_doubles = unique(exc_list_doubles,'rows');
    

end