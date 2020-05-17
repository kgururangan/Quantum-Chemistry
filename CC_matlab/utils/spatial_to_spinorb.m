function [xout] = spatial_to_spinorb(xin)

    Norb = size(xin,1);
    ndim = length(size(xin));
    
    
    if ndim == 2
        
        xout = zeros(2*Norb,2*Norb);
        
        for p = 1:2*Norb
            for q = 1:2*Norb
                if mod(p,2) == mod(q,2)
                    xout(p,q) = xin(floor((p+1)/2),floor((q+1)/2));
                end
            end
        end
        
    else
        
        xout = zeros(2*Norb,2*Norb,2*Norb,2*Norb);
        
        for p = 1:2*Norb
            for q = 1:2*Norb
                for r = 1:2*Norb
                    for s = 1:2*Norb
                        if mod(p,2) == mod(r,2) && mod(q,2) == mod(s,2)
                            xout(p,q,r,s) = xin(floor((p+1)/2),floor((q+1)/2),...
                                                floor((r+1)/2),floor((s+1)/2));
                        end
                    end
                end
            end
        end
        
    end
    
    

end

