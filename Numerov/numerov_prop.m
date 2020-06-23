function [psi] = numerov_prop(x,k2,dir)
    %Nx = floor(abs(Rlim(2)-Rlim(1))/dx)+1;
    Nx = length(x); dx = abs(x(2)-x(1)); dx2 = dx^2;
    psi = zeros(1,Nx);
    %dx = abs(Rlim(2)-Rlim(1))/(Nx-1);
    fn = @(i) 1+dx2/12*k2(i); % FIX k2(i) HERE!!! NOT THE RIGHT ONE FOR LEFT AND RIGHT BOUNDARIES!!!
    if strcmp(dir,'left')
        psi(1) = 0; psi(2) = 1e-20;
        for i = 2:Nx-1
            psi(i+1) = ( (12-10*fn(i))*psi(i) - fn(i-1)*psi(i-1) )/fn(i+1);
        end
    else
        psi(end) = 1e-30; % arbitrary small number > 0 since we are not at +infinity
        psi(end-1) = exp(x(end-1)*sqrt(-k2(end-1)))/exp(x(end)*sqrt(-k2(end)))*psi(end); % JWKB formula
        for i = Nx-1:-1:2
            psi(i-1) = ( (12-10*fn(i))*psi(i) - fn(i+1)*psi(i+1) )/fn(i-1);
        end

    end
end