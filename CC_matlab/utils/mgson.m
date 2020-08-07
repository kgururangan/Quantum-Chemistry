function [Q2, Q, R] = mgson(X, varargin)
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
R = zeros(m,n);
Q = zeros(d,m);
Q2 = [];
for i = 1:m
    v = X(:,i);
    for j = 1:i-1
        R(j,i) = Q(:,j)'*v;
        v = v-R(j,i)*Q(:,j);
    end
    R(i,i) = norm(v);
    Q(:,i) = v/R(i,i);
    if norm(v) > thresh_vec
        Q2(:,i) = v/R(i,i);
    end
end
R(:,m+1:n) = Q'*X(:,m+1:n);

end