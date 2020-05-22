function [VVmat] = eri_fft(orbs,L,h,lambda)

    dx = h(1); dy = h(2); dz = h(3);

    Nx = length(-L(1)/2:dx:L(1)/2);
    Ny = length(-L(2)/2:dy:L(2)/2);
    Nz = length(-L(3)/2:dz:L(3)/2);

    Reff = sqrt(3)*max(L); % needs to be sqrt(3)*L not sqrt(3)*L/2

    % FFTaxis return linear frequency  
    [kx, ~, ~] = FFTaxis( Nx, dx, 1 ); kx = kx*(2*pi); % rad/length
    [ky, ~, ~] = FFTaxis( Ny, dy, 1 ); ky = ky*(2*pi); % rad/length
    [kz, ~, ~] = FFTaxis( Nz, dz, 1 ); kz = kz*(2*pi); % rad/length
    
    % calculate analytical coulomb potential
    Gk = zeros(Nx,Ny,Nz);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz
                kr = sqrt(kx(i)^2+ky(j)^2+kz(k)^2);
                if kr == 0
                    Gk(i,j,k) = Reff^2/2;
                else
                    Gk(i,j,k) = 4*pi*(1-cos(kr*Reff))/(kr^2+lambda^2);
                end
            end
        end
    end
    Gk = fftshift(fftshift(fftshift(Gk,1),2),3); % IMPORTANT: fftshift analytical FT to match with MATLAB's stupid-ass conventions


    Norb = length(orbs);
    VVmat = zeros(Norb,Norb,Norb,Norb);
    
    % loop through permutationally unique 2-body integrals
    for p = 0:Norb-1
        for r = 0:p
            pr = p*(p+1)/2+r;
            for q = 0:Norb-1
                for s = 0:q
                    qs = floor(q*(q+1)/2)+s;
                    if pr >= qs % could we need pr >= qs??
                      
                        fqs = orbs{q+1}.*orbs{s+1};
                        fpr = orbs{p+1}.*orbs{r+1};
                        
                        fqs_k = fftn(fqs); % fourier transform of integral product

                        % calcualte product of fft's equivalent to convolution
                        Gamma_k = fqs_k.*Gk;

                        % ifft the convolution product back to real space
                        Gamma = ifftn(Gamma_k);

                        % calculate integrand as fpr.*ifft(fft(coul).*fft(fqs))
                        integrand = fpr.*Gamma;

                        % calculate integral by simple rectangle rule
                        val = real(sum(sum(sum(integrand,1)*dx,2)*dy,3)*dz);
                        
                        VVmat(p+1,q+1,r+1,s+1) = val;
                        VVmat(q+1,p+1,s+1,r+1) = val;
                        VVmat(r+1,s+1,p+1,q+1) = val;
                        VVmat(s+1,r+1,q+1,p+1) = val;
                        VVmat(r+1,q+1,p+1,s+1) = val;
                        VVmat(p+1,s+1,r+1,q+1) = val;
                        VVmat(s+1,p+1,q+1,r+1) = val;
                        VVmat(q+1,r+1,s+1,p+1) = val;
                        
                    end
                    
                end
            end
        end
    end
    
end
    
    



