function [f] = boys_debug(n,T)
% Exact value of boys function $F_n(T) = \int_{0}^{1} dx x^{2n} exp(-Tx^2)$
% evaluated using Kummer's confluent hypergeometric function 1F1(a,b,x)
    f = hypergeom(n+0.5,n+1.5,-T)/(2*n+1);
end

