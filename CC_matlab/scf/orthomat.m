function [X] = orthomat(S)

    [U, evalS] = eig(S);
    diagS_minushalf = diag(diag(evalS).^(-0.5));
    X = U*diagS_minushalf*U';

end