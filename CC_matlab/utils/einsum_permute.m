function [Aout] = einsum_permute(varargin)

    
    if ischar(varargin{1})
        str = varargin{1};
        Ain = varargin{2};
    else
        str = varargin{2};
        Ain = varargin{1};
    end


    L = split(str,'-');
    IN = strip(L{1});
    OUT = L{2}; OUT = strip(OUT(2:end));

    n = length(IN); m = length(OUT);
    if n ~= m

      error(['cannot permute 1st argument of dimension %d '...
            ' into 2nd argument of dimension %d'],...
            n, m)

    else

        dimA = length(size(Ain)); 
        if n ~= dimA

              error(['contraction string %s does not match rank of input tensor'],str)

        else

            permidx = zeros(1,n);

            for j = 1:n

                idx = find(IN == OUT(j));
                permidx(j) = idx;

            end
    
        end


    end

    Aout = permute(Ain,permidx);

end

