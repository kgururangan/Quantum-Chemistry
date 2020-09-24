function [Q] = mgson_v2(X, varargin)
% Modified Gram-Schmidt orthonormalization (numerical stable version of Gram-Schmidt algorithm) 
% which produces the same result as [Q,R]=qr(X,0)
% Written by Mo Chen (sth4nth@gmail.com).

if isempty(varargin)
    thresh_vec = 1e-15;
else
    thresh_vec = varargin{1};
end

[d,n] = size(X);
m = min(d,n);
Q = [];
Q(:,1) = X(:,1);
for i = 1:m
    v = X(:,i);
    for j = 1:size(Q,2)
        v = v-(Q(:,j)'*v)*Q(:,j);
    end
    if norm(v) > thresh_vec
        Q(:,i) = v;
    end
end

end