function [Hiajb] = HBar_mat_singles(i,a,j,b,HBar)

    Hiajb = HBar{2}{2,1,1,2}(a,j,i,b) + HBar{1}{2,2}(a,b)*deltafcn(i,j) - HBar{1}{1,1}(j,i)*delta(a,b);

end

