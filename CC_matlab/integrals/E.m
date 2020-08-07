function [xout] = E(i,j,t,Rab,a,b)

    p = a + b;
    q = a*b/p;
    if ( t < 0 ) || ( t > (i+j) )
        xout = 0.0;
    elseif i == j && j == t && i == 0
        xout = exp(-q*Rab*Rab);
    elseif j == 0
        xout = (1/(2*p))*E(i-1,j,t-1,Rab,a,b) - (q*Rab/a)*E(i-1,j,t,Rab,a,b) +...
                (t+1)*E(i-1,j,t+1,Rab,a,b);
    else
        xout = (1/(2*p))*E(i,j-1,t-1,Rab,a,b) + (q*Rab/b)*E(i,j-1,t,Rab,a,b) +...
                (t+1)*E(i,j-1,t+1,Rab,a,b);
    end

end

