
M = 16;
L = 1e6;
ns = 10;

z = linspace(0,L,ns*L);

boysval = zeros(1,length(z));
for i = 1:length(z)
    i
    boysval(i) = boys(M,z(i));
end

% n = 1:1:M;
% z = linspace(0,L,ns*L);
% 
% boysval = zeros(length(n),length(z));
% for j = 1:length(n)
%     tic
%     fprintf('\nn = %d',n(j))
%     for k = 1:length(z)
%         boysval(j,k) = boys(n(j),z(k));
%     end
%     fprintf('\nElapsed time - %4.2f s',toc)
% end
% 
% save('boysfunction_values','boysval')