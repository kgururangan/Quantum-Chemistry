function [p,Rp,K] = gauss_prod(alpha,Ra,beta,Rb)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    p = alpha+beta;
    Rp = (alpha*Ra + beta*Rb)/p;
    K = exp(-alpha*beta/p*sum( (Ra-Rb).^2 ));

end

