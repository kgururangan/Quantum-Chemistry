function [g] = gaussian_orbital(g_struct, r, dim)

%     r = linspace(xrange(1),xrange(2),n);
    g = zeros(1,length(r));

    c = g_struct.coeff;
    xi = g_struct.exps;
    R = g_struct.origin;
    R0 = R(dim);

    
    
    nG = length(xi);
    for i = 1:nG
        g = g + c(i)*exp(-xi(i)*(r-R0).^2);
    end

end

