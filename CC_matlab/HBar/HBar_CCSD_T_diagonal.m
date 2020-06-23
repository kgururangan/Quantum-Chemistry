function [D] = HBar_CCSD_T_diagonal(a,b,c,i,j,k,HBar)

%     if length(idx_p) == 3
%         a = idx_p(1); b = idx_p(2); c = idx_p(3);
%         i = idx_h(1); j = idx_h(2); k = idx_h(3);
        D = triples_onebody_diagonal(a,b,c,i,j,k,HBar) + ...
                  triples_twobody_diagonal(a,b,c,i,j,k,HBar) + ...
                  triples_threebody_diagonal(a,b,c,i,j,k,HBar);
% %     elseif length(idx_p) == 2
% %         a = idx_p(1); b = idx_p(2);
% %         i = idx_h(1); j = idx_h(2);
% %         Dabij = doubles_onebody_diagonal(a,b,i,j,HBar) + ...
% %                 doubles_twobody_diagonal(a,b,i,j,HBar);
%     elseif length(idx_p) == 1
%         a = idx_p(1); 
%         i = idx_h(1);
%         D = singles_onebody_diagonal(a,i,HBar) + ...
%               singles_twobody_diagonal(a,i,HBar);
%     end

end

function [Dai] = singles_onebody_diagonal(a,i,HBar)
    Dai = HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i);
end

function [Dai] = singles_twobody_diagonal(a,i,HBar)
    Dai = HBar{2}{2,1,1,2}(a,i,i,a);
end

% function [Dabij] = doubles_onebody_diagonal(a,b,i,j,HBar)
%     Dabij = HBar{1}{2,2}(a,a) + HBar{1}{2,2}(b,b) - ...
%             HBar{1}{1,1}(i,i) - HBar{1}{1,1}(j,j);
% end
%     
% function [Dabij] = doubles_twobody_diagonal(a,b,i,j,HBar)
%     Dabij = HBar{2}{2,1,1,2}(a,i,i,a);
% end
% 
% function [Dabij] = doubles_threebody_diagonal(a,b,i,j,HBar)
% 
% end

function [Dabcijk] = triples_onebody_diagonal(a,b,c,i,j,k,HBar)
    Dabcijk = HBar{1}{2,2}(a,a) + HBar{1}{2,2}(b,b) + HBar{1}{2,2}(c,c) - ...
              HBar{1}{1,1}(i,i) - HBar{1}{1,1}(j,j) - HBar{1}{1,1}(k,k);
end

function [Dabcijk] = triples_twobody_diagonal(a,b,c,i,j,k,HBar)
    Dabcijk = -HBar{2}{2,1,2,1}(a,i,a,i) - HBar{2}{2,1,2,1}(b,i,b,i) - HBar{2}{2,1,2,1}(c,i,c,i)...
              -HBar{2}{2,1,2,1}(a,j,a,j) - HBar{2}{2,1,2,1}(b,j,b,j) - HBar{2}{2,1,2,1}(c,j,c,j)...
              -HBar{2}{2,1,2,1}(a,k,a,k) - HBar{2}{2,1,2,1}(b,k,b,k) - HBar{2}{2,1,2,1}(c,k,c,k)...
              +HBar{2}{1,1,1,1}(i,j,i,j) + HBar{2}{1,1,1,1}(i,k,i,k) + HBar{2}{1,1,1,1}(j,k,j,k)...
              +HBar{2}{2,2,2,2}(a,b,a,b) + HBar{2}{2,2,2,2}(a,c,a,c) + HBar{2}{2,2,2,2}(b,c,b,c);
end

function [Dabcijk] = triples_threebody_diagonal(a,b,c,i,j,k,HBar)
    Dabcijk = -HBar{3}{2,1,1,1,2,1}(a,i,j,i,a,j) - HBar{3}{2,1,1,1,2,1}(a,i,k,i,a,k) - HBar{3}{2,1,1,1,2,1}(a,j,k,j,a,k)...
              -HBar{3}{2,1,1,1,2,1}(b,i,j,i,b,j) - HBar{3}{2,1,1,1,2,1}(b,i,k,i,b,k) - HBar{3}{2,1,1,1,2,1}(b,j,k,j,b,k)...
              -HBar{3}{2,1,1,1,2,1}(c,i,j,i,c,j) - HBar{3}{2,1,1,1,2,1}(c,i,k,i,c,k) - HBar{3}{2,1,1,1,2,1}(c,j,k,j,c,k)...
              +HBar{3}{2,1,2,1,2,2}(a,i,b,i,a,b) + HBar{3}{2,1,2,1,2,2}(a,i,c,i,a,c) + HBar{3}{2,1,2,1,2,2}(b,i,c,i,b,c)...
              +HBar{3}{2,1,2,1,2,2}(a,j,b,j,a,b) + HBar{3}{2,1,2,1,2,2}(a,j,c,j,a,c) + HBar{3}{2,1,2,1,2,2}(b,j,c,j,b,c)...
              +HBar{3}{2,1,2,1,2,2}(a,k,b,k,a,b) + HBar{3}{2,1,2,1,2,2}(a,k,c,k,a,c) + HBar{3}{2,1,2,1,2,2}(b,k,c,k,b,c);
end

