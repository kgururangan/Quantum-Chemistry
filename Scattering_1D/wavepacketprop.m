%****************************************************************
% Program 6: Wavepacket propagation using exponential of H
%****************************************************************
% Parameters for solving the problem in the interval 0 < x < L
L = 100; % Interval Length
N = 400; % No of points
x = linspace(0,L,N)'; % Coordinate vector
dx = x(2) - x(1); % Coordinate step
% Parameters for making intial momentum space wavefunction phi(k)
ko = 2; % Peak momentum
a = 20; % Momentum width parameter
dk = 2*pi/L; % Momentum step
km=N*dk; % Momentum limit
k=linspace(0,+km,N)'; % Momentum vector
% Make psi(x,0) from Gaussian kspace wavefunction phi(k) using
% fast fourier transform :
phi = exp(-a*(k-ko).^2).*exp(-1i*6*k.^2); % unnormalized phi(k)
psi = ifft(phi); % multiplies phi by expikx and integrates vs. x
psi = psi/sqrt(psi'*psi*dx); % normalize the psi(x,0)
% Expectation value of energy; e.g. for the parameters
% chosen above <E> = 2.062.
avgE = phi'*0.5*diag(k.^2,0)*phi*dk/(phi'*phi*dk);
% CHOOSE POTENTIAL U(X): Either U = 0 OR
% U = step potential of height avgE that is located at x=L/2
%U = 0*heaviside(x-(L/2)); % free particle wave packet evolution
U = avgE*heaviside(x-(L/2)); % scattering off step potential
% Finite-difference representation of Laplacian and Hamiltonian
e = ones(N,1); Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
H = -(1/2)*Lap + spdiags(U,0,N,N);
% Parameters for computing psi(x,T) at different times 0 < T < TF
NT = 200; % No. of time steps
TF = 29; T = linspace(0,TF,NT); % Time vector
dT = T(2)-T(1); % Time step
hbar = 1;
% Time displacement operator E=exp(-iHdT/hbar)
E = expm(-1i*full(H)*dT/hbar); % time desplacement operator
%***************************************************************
% Simulate rho(x,T) and plot for each T
%***************************************************************
for t = 1:NT; % time index for loop
% calculate probability density rho(x,T)
psi = E*psi; % calculate new psi from old psi
rho = conj(psi).*psi; % rho(x,T)
plot(x,rho,'k'); % plot rho(x,T) vs. x
axis([0 L 0 0.15]); % set x,y axis parameters for plotting
xlabel('x [m]', 'FontSize', 24);
ylabel('probability density [1/m]','FontSize', 24);
pause(0.05); % pause between each frame displayed
end
% Calculate Reflection probability
R=0;
for a=1:N/2
R=R+rho(a);
end
R=R*dx;