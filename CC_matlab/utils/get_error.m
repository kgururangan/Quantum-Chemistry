function [err] = get_error(A,B)
    % quick routine to get error between two tensors
    err = A - B;
    err = sum(abs(err(:)));
end

