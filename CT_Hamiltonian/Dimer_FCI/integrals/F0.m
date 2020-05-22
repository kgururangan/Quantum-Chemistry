function [xout] = F0(t)
    if sqrt(t) <= eps
        xout = 1 - t/3;
    else
        xout = 1/2*sqrt(pi/t)*erf(sqrt(t));
    end
end