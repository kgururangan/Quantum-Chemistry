function [xout] = overlap(a,lmn1,A,b,lmn2,B)
    Sx = E(lmn1(1),lmn2(1),0,A(1)-B(1),a,b); % X
    Sy = E(lmn1(2),lmn2(2),0,A(2)-B(2),a,b); % Y
    Sz = E(lmn1(3),lmn2(3),0,A(3)-B(3),a,b); % Z
    xout = Sx*Sy*Sz*power(pi/(a+b),1.5);
end

