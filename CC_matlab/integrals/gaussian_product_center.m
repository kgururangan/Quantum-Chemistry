function [p,P] = gaussian_product_center(a,A,b,B)
    p = a + b;
    P = (a*A + b*B)/p;
end

