function [L, R] = biorth(L, R)

% biorthogonalizes two matrices L and R such that L*R = I
% note dimensions should be such that dim(R') = dim(L)

% M = L*R 
% L*R = M -> L*R = ML*MU using LU decomposition
% invert since L and U matrices are simpler to invert
% inv(ML)*L*R*inv(MU) = I = Lp*Rp
% therefore, 
% Lp = inv(ML)*L -> Lp = ML\L
% Rp = R*inv(MU) -> Rp = R/MU

%     if size(R') ~= size(L)
%         L = L';
%     end

    LR = L*R;
    if isdiag(LR)
        d = diag(LR);
        Msqrt = diag(sqrt(d));
        L = L/Msqrt;
        R = R/Msqrt;
    else
        [ML,MU] = lu(LR);
         L = ML\L;
         R = R/MU;
    end
    
end

% function [ V, W ] = biorth( V, W )
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% 
%     W = mgson(W);
% 
%     for i = 1:size(V,2)
%         v = V(:,i);
%         for j = 1:size(W,2)
%             if i ~= j
%                 w = W(:,j)/norm(W(:,j));
%                 v = v - (ctranspose(w)*v)*w;
%             end
%         end
%         v = v/(ctranspose(W(:,i))*v);
%         V(:,i) = v;
%     end
%         
% 
% 
% end

