
clear all
%clc
close all
%

h = 0.1;

dx = h; dy = h; dz = h;

L = 10;

fprintf('WRAPPER...\n')
xi1 = 1.0; xi2 = 1.0; xi3 = 2.0; xi4 = 2.0;
R1 = [-1,0,0]; R2 = [2,0,0];

Rat(1,:) = R1; Rat(2,:) = R2;
g1.origin = R1; g2.origin = R2; g3.origin = R1; g4.origin = R2;
g1.exps = xi1; g2.exps = [xi2,xi4]; g3.exps = xi3; g4.exps = xi4;
g1.coeff = 1.0; g2.coeff = [0.5, 0.5]; g3.coeff = 1.0; g4.coeff = 1.0;
orbs = {g1,g2,g3,g4};

[error, VVmat_fft, VVmat] = test_eri_fft(orbs,Rat,[L,L,L],[h,h,h]);

%%
xi1 = 1.0; xi2 = 1.0; xi3 = 2.0; xi4 = 2.0;
R1 = [-1,0,0]; R2 = [1,0,0];

Rat(1,:) = R1; Rat(2,:) = R2;
g1.origin = R1; g2.origin = R2; g3.origin = R1; g4.origin = R2;
g1.exps = [xi1,xi2,0]; g2.exps = [xi2,xi3,0]; g3.exps = [xi1,xi3,0]; g4.exps = [xi2,xi4,0];
g1.coeff = [0.5,0.5,0]; g2.coeff = [0.5,0.5,0]; g3.coeff = [0.5,0.5,0]; g4.coeff = [0.5,0.5,0];
g1.shell = [0,0,0]; g2.shell = [0,0,0]; g3.shell = [0,0,0]; g4.shell = [0,0,0];

[Smat, ~, ~, VVmat] = spatial_integrals_v2({g1,g2,g3,g4},Rat,[1,1]);

S1 = zeros(Nx,Ny,Nz);
S2 = zeros(Nx,Ny,Nz);
S3 = zeros(Nx,Ny,Nz);
S4 = zeros(Nx,Ny,Nz);

for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            
            sh1 = g1.shell; sh2 = g2.shell; sh3 = g3.shell; sh4 = g4.shell;
            
            num_gauss = length(g1.exps);
            
            for I = 1:num_gauss
                
                temp1 =     g1.coeff(I)*(x(i)-g1.origin(1))^sh1(1)*(y(j)-g1.origin(2))^sh1(2)*(z(k)-g1.origin(3))^sh1(3).*...
                            exp(-g1.exps(I)*((x(i)-g1.origin(1))^2+(y(j)-g1.origin(2))^2+(z(k)-g1.origin(3))^2));
                temp2 =     g2.coeff(I)*(x(i)-g2.origin(1))^sh2(1)*(y(j)-g2.origin(2))^sh2(2)*(z(k)-g2.origin(3))^sh2(3).*...
                            exp(-g2.exps(I)*((x(i)-g2.origin(1))^2+(y(j)-g2.origin(2))^2+(z(k)-g2.origin(3))^2));
                temp3 =     g3.coeff(I)*(x(i)-g3.origin(1))^sh3(1)*(y(j)-g3.origin(2))^sh3(2)*(z(k)-g3.origin(3))^sh3(3).*...
                            exp(-g3.exps(I)*((x(i)-g3.origin(1))^2+(y(j)-g3.origin(2))^2+(z(k)-g3.origin(3))^2));
                temp4 =     g4.coeff(I)*(x(i)-g4.origin(1))^sh4(1)*(y(j)-g4.origin(2))^sh4(2)*(z(k)-g4.origin(3))^sh4(3).*...
                            exp(-g4.exps(I)*((x(i)-g4.origin(1))^2+(y(j)-g4.origin(2))^2+(z(k)-g4.origin(3))^2));
                        
                S1(i,j,k) = S1(i,j,k) + temp1*primitive_norm(g1.shell(1),g1.shell(2),g1.shell(3),g1.exps(I));
                S2(i,j,k) = S2(i,j,k) + temp2*primitive_norm(g2.shell(1),g2.shell(2),g2.shell(3),g2.exps(I));
                S3(i,j,k) = S3(i,j,k) + temp3*primitive_norm(g3.shell(1),g3.shell(2),g3.shell(3),g3.exps(I));
                S4(i,j,k) = S4(i,j,k) + temp4*primitive_norm(g4.shell(1),g4.shell(2),g4.shell(3),g4.exps(I));
                
            end
            
        end
    end
end
fprintf('Spatial orbitals computed on mesh\n')

ORBS = {S1,S2,S3,S4};

Norb = length(ORBS);

Smat2 = zeros(Norb);
for i = 1:Norb
    for j = 1:Norb
        Smat2(i,j) = sum(sum(sum(ORBS{i}.*ORBS{j},1)*dx,2)*dy,3)*dz;
    end
end
Smat2

% [X,Y,Z] = meshgrid(x,y,z);
% 
% S1m = exp(-xi1*((X-R1(1)).^2+(Y-R1(2)).^2+(Z-R1(3)).^2))*primitive_norm(0,0,0,xi1);

%%
%
fprintf('Explicit\n')

x = -L/2:dx:L/2;
y = -L/2:dy:L/2;
z = -L/2:dz:L/2;

Reff = sqrt(3)*L;

Nx = length(x); Ny = length(y); Nz = length(z);

% FFTaxis return linear frequency  
[kx, ~, ~] = FFTaxis( length(x), dx, 1 ); %kx = ifftshift(kx)*(2*pi); % rad/length
[ky, ~, ~] = FFTaxis( length(y), dy, 1 ); %ky = ifftshift(ky)*(2*pi); % rad/length
[kz, ~, ~] = FFTaxis( length(z), dz, 1 ); %kz = ifftshift(kz)*(2*pi); % rad/length

kx = kx*2*pi;
ky = ky*2*pi;
kz = kz*2*pi;


for p = 1:2
    for q = 1:2
        for r = 1:2
            for s = 1:2

             [val] = eri_fft_integrate(ORBS{p}, ORBS{q}, ORBS{r}, ORBS{s}, L, h, 0);
             err = abs( (VVmat(p,q,r,s)-val)/VVmat(p,q,r,s) )*100;
             fprintf('Gaussian integration gives V(%d,%d,%d,%d) = %4.6f\n',p,q,r,s,VVmat(p,q,r,s))
             fprintf('FFT integration gives V(%d,%d,%d,%d) = %4.6f\n',p,q,r,s,val)
             fprintf('Error = %4.6f%% \n\n',err)
            end
        end
    end
end
%%


% (S1,S2|S1,S2)
for p = 1:2
    for q = 1:2
        for r = 1:2
            for s = 1:2


                fprintf('Gaussian integration gives V(%d,%d,%d,%d) = %4.6f\n',p,q,r,s,VVmat(p,q,r,s))

                tic
                % calculate fourier transform of coulomb kernel analytically
                % factors of 2*pi.>?
                Gk = zeros(Nx,Ny,Nz);
                G = zeros(Nx,Ny,Nz);
                fpr = zeros(Nx,Ny,Nz);
                fqs = zeros(Nx,Ny,Nz);

                lambda = 0; 
                for i = 1:Nx
                    for j = 1:Ny
                        for k = 1:Nz

                            fqs(i,j,k) = ORBS{q}(i,j,k)*ORBS{s}(i,j,k);
                            fpr(i,j,k) = ORBS{p}(i,j,k)*ORBS{r}(i,j,k);

                            kr = sqrt(kx(i)^2+ky(j)^2+kz(k)^2);
                            if kr == 0
                                Gk(i,j,k) = Reff^2/2;
                                %Gk(i,j,k) = 1;
                            else
                                Gk(i,j,k) = 4*pi*(1-cos(kr*Reff))/(kr^2+lambda^2);
                                %Gk(i,j,k) = 2*pi*(1-cos(kr*Reff))/(kr^2+lambda^2);
                                %Gk(i,j,k) = 1/(kr^2+lambda^2);
                            end



                        end
                    end
                end

                %
                Gk = fftshift(fftshift(fftshift(Gk,1),2),3);
                fqs_k = fftn(fqs); % fourier transform of integral product

                % calcualte product of fft's equivalent to convolution
                Gamma_k = zeros(Nx,Ny,Nz);
                for i = 1:Nx
                    for j = 1:Ny
                        for k = 1:Nz
                            Gamma_k(i,j,k) = fqs_k(i,j,k)*Gk(i,j,k);
                        end
                    end
                end

                % ifft the convolution product back to real space
                Gamma = ifftn(Gamma_k);%/(dx*dy*dz); 


                % calculate integrand as fpr.*ifft(fft(coul).*fft(fqs))
                temp = zeros(Nx,Ny,Nz);

                %
                for i = 1:Nx
                    for j = 1:Ny
                        for k = 1:Nz
                            temp(i,j,k) = fpr(i,j,k)*Gamma(i,j,k);
                        end
                    end
                end
                figure(3)
                imagesc(x,y,squeeze(real(temp(:,:,floor(Nz/2)))))
                axis square
                colorbar
                title('Integrand')

                % calculate integral by simple rectangle rule
                INT = sum(sum(sum(temp,1)*dx,2)*dy,3)*dz;
                %INT = sum(sum(sum(temp,1)*dx,2)*dy,3)*dz;

                fprintf('FFT integration gives V(%d,%d,%d,%d) = %4.6f (%4.1f seconds)\n',p,q,r,s,INT,toc)

            end
        end
    end
end

%%

subplot(211)
imagesc(x,y,squeeze(fqs(:,:,floor(Nz/2))))
axis square
colorbar

fqs_kk = fft_shift(fqs);

subplot(212)
imagesc(kx,ky,squeeze(abs(fqs_kk(:,:,floor(Nz/2)))))
axis square
colorbar




%%
[val] = eri_fft_integrate(S1, S1, S2, S2, L, h, 0)

%%

for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            TEMP(i,j,k) = fpr(i,j,k)*G(i,j,k)*fqs(i,j,k);
        end
    end
end
INT = sum(sum(sum(TEMP,1)*dx,2)*dy,3)*dz

%%

lambda = 0.01;
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            Gk(i,j,k) = 1/( (kx(i)^2 + ky(j)^2 + kz(k)^2) + lambda^2 );
            if x(i)^2 + y(j)^2 + z(k)^2 ~= 0
                G(i,j,k) = 1/(sqrt(x(i)^2+y(j)^2+z(k)^2));
            else
                G(i,j,k) = -lambda;
            end
        end
    end
end

Gk2 = fft_shift(G);



%%
clf
figure(1)
for mm = 1:length(RJ)
    
subplot(211)
imagesc(XI1,XI2,squeeze(Int_fft(:,:,mm)))
axis square
colorbar

subplot(212)
imagesc(XI1,XI2,squeeze(Exact(:,:,mm)))
axis square
colorbar

pause(2)

end
%%
A = rand(10); B = rand(10);
CONV = reshape(fftshift(cconv(fqs(:),G(:),Nx*Ny*Nz)),Nx,Ny,Nz);

figure(4)
imagesc(x,y,squeeze(real(CONV(:,:,floor(Nz/2)))))
axis square
colorbar
title('Integrand')
%%

a = exp(-(x-1).^2);
b = exp(-(x-0.5).^2);

c1 = cconv(a,b,Nx);
c2 = fftshift(fft(a)).*fftshift(fft(b));
c2 = ifft(c2,Nx);

plot(x,c1,x,c2)

%% Analytical stuff

xi1 = 1; R1 = [1,0,0];
xi2 = 1; R2 = [-1,0,0];
[p,Rp,C] = gauss_prod(xi1,R1,xi2,R2);

analytical = zeros(Nx,Ny,Nz);
S = zeros(Nx,Ny,Nz);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            analytical(i,j,k) = primitive_norm(0,0,0,xi1)*primitive_norm(0,0,0,xi2)*C*...
                         exp(-1/(4*p)*(kx(i)^2+ky(j)^2+kz(k)^2) - 1i*(kx(i)*Rp(1)+ky(j)*Rp(2)+kz(k)*Rp(3)));
            S1(i,j,k) = primitive_norm(0,0,0,xi1)*primitive_norm(0,0,0,xi2)*C*...
                        exp(-p*((x(i)-Rp(1))^2+(y(j)-Rp(2))^2+(z(k)-Rp(3))^2));
        end
    end
end

num = fftn(S1);
S2 = ifftn(num);

numerical = fft_shift(S1)*dx*dy*dz;


%%
figure(4)
subplot(211)
imagesc(kx,ky,squeeze(real(S1(:,:,floor(Nz/2)))))
axis square
colorbar
subplot(212)
imagesc(kx,ky,squeeze(real(S2(:,:,floor(Nz/2)))))
axis square
colorbar
%%
clear all

dx = 0.01;
x = -3:dx:3;

kx = FFTaxis_kg(x,1,0,0);

alpha = 10;
x0 = 1;

y = exp(-alpha*(x-x0).^2)*primitive_norm(0,0,0,alpha);

yf = compute_dft_vectorized(y);
y1 = (fft(y));

plot(kx,yf,kx,y1)

% y1 = fftshift(fft(y))*dx;%*exp(-alpha*x0^2);
% 
% [frq,y1R,y1I,ampy1,~] = simpleFFT(y,dx);
% 
% y2 = sqrt(pi/alpha)*exp(-alpha*x0^2).*exp(-1/(4*alpha)*(kx+2*1i*alpha*x0).^2);
% 
% plot(kx,y1,kx,y2)


%%
fs = 500;
dt = 1/fs;
tmax = 4;
t = 0:dt:tmax;

f1 = [10, 20, 30, 40, 50, 60];
c = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];

signal = zeros(1,length(t));
for i = 1:length(f1)
    signal = signal + c(i)*sin(2*pi*f1(i)*t);
end

[frq, real_part, imag_part, amp, phase] = simpleFFT( signal, fs);

plot(frq,amp)
%%
dx = 0.001;
x = -10:dx:10;
Nx = length(x);

a = 1;
A = 1;
rect1 = zeros(1,Nx); 
rect2 = zeros(1,Nx);
tri1 = zeros(1,Nx);

for i = 1:Nx
    if abs(x(i)) <= a
        rect1(i) = A;
        rect2(i) = A;
    end
    if abs(x(i)) <= 2*a
        if x(i) <= 0
            tri1(i) = A^2*(x(i)+2*a);
        else
            tri1(i) = A^2*(2*a-x(i));
        end
    end
end

tri2 = fftshift(cconv(rect1,rect2,Nx))*dx;
tri3 = ifftshift(ifft(fft(rect1).*fft(rect2)))*dx;

plot(x,tri2,x,tri3)


%%

function output = compute_dft_vectorized(input)
	assert(isvector(input));
	n = numel(input);
	matrix = exp(-2i * pi / n * (0 : n-1)' * (0 : n-1));
	if size(input, 1) == 1  % Row vector
		output = input * matrix;
	elseif size(input, 2) == 1  % Column vector
		output = matrix * input;
	end
end

function output = compute_dft_scalarized(input)
	assert(isvector(input));
	n = numel(input);
	output = NaN(size(input));
	for k = 0 : n - 1  % For each output element
		s = 0;
		for t = 0 : n - 1  % For each input element
			s = s + input(t + 1) * exp(-2i * pi * t * k / n);
		end
		output(k + 1) = s;
	end
end


function y = fft_shift(x)
    y = fftshift(fftshift(fftshift(fftn(x),1),2),3);
end


function y = ifft_shift(x)
    y = fftshift(fftshift(fftshift(ifftn(x),1),2),3);
end
