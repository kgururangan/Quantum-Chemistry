function [sign] =  signPermutation(p)
% ----------------------------------------------------------
% Calculates the sign of a permutation p.
% p is a row vector p(1,n), which represents the permutation.
% sgn(p) = (-1)^(No. of even-length cycles)
% Complexity : O(n + ncyc) ~ O(n + Hn) ~~ O(n+log(n)) steps.
%
% Derek O'Connor 20 March 2011.
% ----------------------------------------------------------

    n   = length(p);
    visited(1:n) = false;                  % Logical vector which marks all p(k)
                                           % not visited
    sign = 1;
    for k = 1:n
        if ~visited(k)                     % k not visited, start of new cycle
            ct = k;
            L = 0;
            while ~visited(ct)             % Traverse the current cycle k
                L = L+1;                   % and find its length L
                visited(ct) =  true;
                ct    = p(ct);
            end
            if mod(L,2) == 0               % If L is even, change sign.
                sign = -1*sign;
            end
        end % if ~visited(k)
    end % for k

end
