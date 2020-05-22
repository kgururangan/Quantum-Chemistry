function [taus, tau] = time_from_freq(omega, c)

    dom = abs(omega(2)-omega(1));
    Nw = length(omega);
    
    tau = 2*pi*(0:1/2/dom/(Nw):(Nw-1)*1/2/dom/Nw); % [cm]
    
    taus = tau/c;  % [fs]
    
end
    

