function [D1,D2,D3] = HBar_CCSD_D_diagonal(a,b,i,j,HBar)

    D1 = HBar{1}{2,2}(a,a) + HBar{1}{2,2}(b,b) - HBar{1}{1,1}(i,i) - HBar{1}{1,1}(j,j);
    D2 = HBar{2}{2,2,2,2}(a,b,a,b) + HBar{2}{1,1,1,1}(i,j,i,j)...
         + HBar{2}{2,1,1,2}(a,i,i,a) + HBar{2}{2,1,1,2}(b,j,j,b) ...
         - HBar{2}{2,1,1,2}(b,i,i,b) - HBar{2}{2,1,1,2}(a,j,j,a);
    D3 = 0;

end
