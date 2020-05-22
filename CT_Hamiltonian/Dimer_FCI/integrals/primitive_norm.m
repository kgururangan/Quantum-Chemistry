function [N] = primitive_norm(l,m,n,alpha)

    N = (2/pi)^(3/4) * (2^(l+m+n)*alpha^((2*l+2*m+2*n+3)/4))/...
                       (sqrt(fact2(2*l-1)*fact2(2*m-1)*fact2(2*n-1)));

end

function xout = fact2(n)

    if n == 0 || n == -1
        xout = 1;
    elseif mod(n,2) == 1
        xout = prod([1:2:n]);
    else
        xout = prod([2:2:n]);
    end
    
end
    
    