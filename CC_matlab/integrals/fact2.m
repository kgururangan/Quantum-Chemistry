function [xout] = fact2(n)
    if mod(n,2) == 0
        xout = prod(2:2:n);
    else
        xout = prod(1:2:n);
    end
end


