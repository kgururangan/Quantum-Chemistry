function [xout] = vABCD(alpha,rA,beta,rB,gamma,rC,delta,rD)
% THESE ARE IN CHEMIST' NOTATION!!!
    p = alpha + beta;
    q = gamma + delta;
    
    rP = (alpha*rA + beta*rB)/p;
    rQ = (gamma*rC + delta*rD)/q;
    
    x1 = alpha*beta/p;
    x2 = gamma*delta/q;
    x3 = sum( (rA - rB).^2 );
    x4 = sum( (rC - rD).^2 );
    x5 = sum( (rP - rQ).^2 );
    x6 = p*q/(p+q)*x5;
    
    xout = 2*pi^(5/2)/(p*q*sqrt(p+q))*exp(-x1*x3 - x2*x4)*F0(x6);

    xout = (2*alpha/pi)^(3/4)*(2*beta/pi)^(3/4)*(2*gamma/pi)^(3/4)*(2*delta/pi)^(3/4)*xout;

end