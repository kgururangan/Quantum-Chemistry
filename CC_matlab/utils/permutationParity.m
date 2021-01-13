function [p] = permutationParity(P,Dim)

[nRows,nCols] = size(P);

if nargin < 2 || isempty(Dim)
    if nRows == 1
        Dim = 2;
    elseif nCols == 1
        Dim = 1;
    else
        Dim = 1;
    end
end

p = 0;

if Dim == 1
    for i = 1:nRows
        p = sum(repmat(P(i,:),nRows-i,1) > P((i+1):end,:),1) + p;
    end
elseif Dim == 2
    for i = 1:nCols
        p = sum(repmat(P(:,i),1,nCols-i) > P(:,(i+1):end),2) + p;
    end
end
p = mod(p,2);
if p == 0
    p = 1;
else
    p = -1;
end
end