function [ V, W ] = biorth( V, W )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    W = mgson(W);

    for i = 1:size(V,2)
        v = V(:,i);
        for j = 1:size(W,2)
            if i ~= j
                w = W(:,j)/norm(W(:,j));
                v = v - (ctranspose(w)*v)*w;
            end
        end
        v = v/(ctranspose(W(:,i))*v);
        V(:,i) = v;
    end
        


end

