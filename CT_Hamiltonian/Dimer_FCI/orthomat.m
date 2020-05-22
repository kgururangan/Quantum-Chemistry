function [X] = orthomat(S,kind)

    [U, evalS] = eig(S);
    
    diagS_minushalf = diag((diag(evalS).^(-0.5)));
    
    if strcmp(kind,'symmetric')
        X = U*(diagS_minushalf)*U';
    else
        X = U*diagS_minushalf;
    end
    
end
