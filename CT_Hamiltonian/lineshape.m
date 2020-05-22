function [I] = lineshape(omega, wab, gamma, type)


if strcmp(type,'lorentzian')
    %I = 1i./(omega - wab - 1i*gamma);
    %I = 1./(gamma - 1i*(omega - (wab)));
    I = 1./(gamma - 1i*(omega-wab));
    %I = I./trapz(omega,I);
    %I = I./norm(I);
end

if strcmp(type,'delta')
    [~,idx] = min(abs(omega-wab));
    I = zeros(length(omega));
    I(idx) = 1;
end

if strcmp(type,'gaussian')
  %  I = exp(-1i*(omega - wab)).*exp(-(omega - wab).^2./(2*gamma^2));
  % I = 1/(gamma*sqrt(2*pi))*exp(-(omega-wab).^2/(2*gamma^2));
    I = exp(-(omega-wab).^2/(gamma^2));
end

if strcmp(type,'dampsine')
    I = exp(1i*wab*omega - gamma*omega);
end

if strcmp(type,'voigt')
    hom_inhom = 1;
    I = voigt(omega,wab,gamma*hom_inhom, gamma);
end
    

end