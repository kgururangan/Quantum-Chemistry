function [xout] = tAB(alpha,rA,beta,rB)
    x1 = alpha*beta/(alpha+beta);
    x2 = sum( (rA-rB).^2 );
    xout = x1*(3-2*x1*x2)*(pi/(alpha+beta))^(3/2)*exp(-x1*x2);
    xout = (2*alpha/pi)^(3/4)*(2*beta/pi)^(3/4)*xout;
end