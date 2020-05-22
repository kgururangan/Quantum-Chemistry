function [norm_coef] = normalize_AO(g_struct)

    xi = g_struct.exps;
    cm = g_struct.coeff;
    origin = g_struct.origin;
    nG = length(xi);
    
    GG = 0.0;
    for i = 1:nG
        for j = 1:nG
            GG = GG + cm(i)*cm(j)*sAB(xi(i),origin,xi(j),origin);
        end
    end
    
    norm_coef  = sqrt(1/GG);
    
  

end

