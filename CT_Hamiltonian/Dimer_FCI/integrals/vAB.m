function [xout] = vAB(alpha,rA,beta,rB,zC,rC)
    p = alpha+beta;
    rP = (alpha*rA + beta*rB)/p;
    x1 = alpha*beta/p;
    x2 = sum( (rA-rB).^2 );
    x3 = sum( (rP-rC).^2 );
    x4 = p*x3;
    
    xout = -2*pi/p*zC*exp(-x1*x2)*F0(x4);
    xout = (2*alpha/pi)^(3/4)*(2*beta/pi)^(3/4)*xout;
end