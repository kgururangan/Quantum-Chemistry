function [val] = eri_fft_integrate(orb_p, orb_q, orb_r, orb_s, L, h, lambda)

% there's a problem here with floating factors of 2.. ad hoc factor of *2
% added to final integration.

% also results don't change much with mesh size (h), but they change
% significantly with box size (L)

    dx = h; dy = h; dz = h;

    x = -L/2:h:L/2;
    y = -L/2:h:L/2;
    z = -L/2:h:L/2;

    Reff = sqrt(3)*L; % needs to be sqrt(3)*L not sqrt(3)*L/2

    Nx = length(x); Ny = length(y); Nz = length(z);

    % FFTaxis return linear frequency  
    [kx, ~, ~] = FFTaxis( length(x), dx, 1 ); kx = kx*(2*pi); % rad/length
    [ky, ~, ~] = FFTaxis( length(y), dy, 1 ); ky = ky*(2*pi); % rad/length
    [kz, ~, ~] = FFTaxis( length(z), dz, 1 ); kz = kz*(2*pi); % rad/length
    
    % NOTE:
    % A centered real space/time axis produces a linear frequency axis
    % (-Nyquist : Nyquist)
    % A time/space axis starting from 0 product a shifted frequency axis
    % (0 : Nyquist, -Nyquist:0)


    % calculate fourier transform of coulomb kernel analytically
    % factors of 2*pi.>?
    Gk = zeros(Nx,Ny,Nz);
    fpr = zeros(Nx,Ny,Nz);
    fqs = zeros(Nx,Ny,Nz);

    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz

                fqs(i,j,k) = orb_q(i,j,k)*orb_s(i,j,k);
                fpr(i,j,k) = orb_p(i,j,k)*orb_r(i,j,k);

                kr = sqrt(kx(i)^2+ky(j)^2+kz(k)^2);
                if kr == 0
                    Gk(i,j,k) = Reff^2/2;
                else
                    Gk(i,j,k) = 4*pi*(1-cos(kr*Reff))/(kr^2+lambda^2);
                end



            end
        end
    end

    %
    
%     fqs = orb_q.*orb_s;
%     fpr = orb_p.*orb_r;
    
    

    fqs_k = fftn(fqs); % fourier transform of integral product
    Gk = fftshift(fftshift(fftshift(Gk,1),2),3); % IMPORTANT: fftshift analytical FT to match with MATLAB's stupid-ass conventions
    
    % calcualte product of fft's equivalent to convolution
    Gamma_k = fqs_k.*Gk;

    % ifft the convolution product back to real space
    Gamma = ifftn(Gamma_k);

    % calculate integrand as fpr.*ifft(fft(coul).*fft(fqs))
    integrand = fpr.*Gamma;

    % calculate integral by simple rectangle rule
    val = real(sum(sum(sum(integrand,1)*dx,2)*dy,3)*dz);

end

