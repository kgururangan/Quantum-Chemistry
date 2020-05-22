function xout = sAB(alpha,rA,beta,rB)
    p = alpha + beta; %x1 = sum( (rA-rB).^2 );
    x1 = dot(rA-rB,rA-rB);
    xout = (pi/p)^(3/2)*exp(-alpha*beta/p*x1);
    xout = (2*alpha/pi)^(3/4)*(2*beta/pi)^(3/4)*xout;
end