function [r] = update_R(r, e, D)

    for i = 1:length(r)
        %fprintf('e = %4.10f\n',e)
        %fprintf('D = %4.10f\n',D(i))
        %fprintf('DELTA = %4.10f\n',e - D(i))
        r(i) = r(i)/(e - D(i));
    end
    
end

