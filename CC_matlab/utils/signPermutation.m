function [sign] =  signPermutation(p)
%  # explicitly using cycle notation of permutation
%  # input: p is a list of permutation indices
%  # i.e. [x_sorted, p] = sort(x_unsorted)

    n = length(p);
    visited = zeros(1,n);
    sign = 1;
    
    for k = 1:n
        if visited(k) == 0
            ct = k;
            L = 0;

            while visited(k) == 0
                    L = L + 1;
                    visited(ct) = 1;
                    ct = p(ct);
            end
            
            if mod(L,2) == 0
                sign = -1*sign;
            end

        end
    end
    
end