function [D1,D2,D3] = HBar_CCSD_S_diagonal(a,i,HBar)

    D1 = HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i);
    D2 = HBar{2}{2,1,1,2}(a,i,i,a);
    D3 = 0;
end

