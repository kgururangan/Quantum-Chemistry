function [ITP,bound_flag] = get_turning_points(V,E)
    [~,ieq] = min(V);
    VL = V(1:ieq); VR = V(ieq:end);
    [~,ITP(1)] = min(abs(VL-E)); [~,ITP(2)] = min(abs(VR-E));
    ITP(2) = ITP(2) + ieq-1;
    if ITP(2) == length(V)
        bound_flag = 1;
    else
        bound_flag = 0;
    end
        
end

