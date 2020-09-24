function [B] = get_diagonal_guess(B,D,nroot)

    I = eye(size(B,1));
    [~,idx] = sort(D,'ascend');
    B(:,1:nroot) = I(:,idx(1:nroot));

end

