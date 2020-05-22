clearvars
clc

% crystal organization: the crystal is assumed to be 1-dimensional with
% alternate donor-acceptor stacking.The odd sites are donors and even sites
% are acceptors so the crystal looks like: DADADADA. Hamiltonian matrix
% organization: the first columns are the single-particle FE states ordered
% by number, i.e. |1> |2> |3> ... where |n> describes a state where there's
% a FE on the n'th crystal site. Then come all the CT states with hole on
% the leftmost donor and going up with increasing electron position on the
% acceptor, i.e. for 2 pairs: |1+,2-> |1+,4-> |3+,2-> |3+,4->.

% all energies are in eV
N_pairs     =   2; % number of donor-acceptor pairs in crystal
ed          =   2.7; % excitation energy of single donor molecule
ea          =   2.25; % excitation energy of single acceptor molecule
vda         =   0.15; % donor-acceptor FE coupling
ectinf      =   1.7; % energy of unbound electron-hole pair (relative to ground state energy)
vctnn       =   1.65; % energy of electron-hole pair residing on nearest-neighbors
te          =   0.05; % electron transfer/dissociation integral
th          =   0.05; % hole transfer/dissociation integral
gamma       =   0.05; % spectral linewidth

% number of states
num_states = 2*N_pairs + N_pairs^2;

% energy axis (eV)
omega = linspace(1, 3.4, 1000);

%%%%%%%%%%%%%%%
% Hamiltonian %
%%%%%%%%%%%%%%%

% H1 = donor FE energy
H1 = blkdiag(diag(repmat([ed,0],1,N_pairs)), zeros(N_pairs^2));

% H2 = acceptor FE energy
H2 = blkdiag(diag(repmat([0,ea],1,N_pairs)), zeros(N_pairs^2));

% H3 = donor-acceptor FE interaction
H3 = blkdiag(diag(vda*ones(1,2*N_pairs-1),1) + diag(vda*ones(1,2*N_pairs-1),-1), zeros(N_pairs^2));

% H4 = charge transfer state energy
ect = @(m) ectinf-1./m*vctnn;
temp = zeros(N_pairs);
for i=1:N_pairs
    for j=1:N_pairs
        temp(i,j) = 2*abs(i-j)+1;
    end
end
temp = ect(temp);
temp = reshape(temp, 1, []);
H4 = blkdiag(zeros(2*N_pairs), diag(temp));
clear temp

% H5 = electron migration energy
temp = diag(te*ones(1,N_pairs-1),1) + diag(te*ones(1,N_pairs-1),-1);
temp(1,N_pairs) = te;
temp(N_pairs,1) = te;
H5 = temp;
for i=1:N_pairs-1
    H5 = blkdiag(H5, temp);
end
H5 = blkdiag(zeros(2*N_pairs), H5);
clear temp

% H6 = hole migration energy
temp = diag(th*ones(1,N_pairs-1),1) + diag(th*ones(1,N_pairs-1),-1);
temp(1,N_pairs) = th;
temp(N_pairs,1) = th;
H6 = kron(temp, eye(N_pairs));
H6 = blkdiag(zeros(2*N_pairs), H6);
clear temp

% H7 = energy for donor-FE dissociation to CT
H7 = zeros(num_states);
temp = [te; zeros(N_pairs-2, 1); te];
for i=1:2:2*N_pairs-1
    j = (i+1)/2;
    H7((2*N_pairs+1)+(j-1)*numel(temp):(2*N_pairs+1)+j*numel(temp)-1, i) = temp;
    temp = circshift(temp,1);
end
H7(2*N_pairs+numel(temp), 1) = 0; % transfer from closed to open boundary conditions
H7 = H7+H7';
clear temp

% H8 = energy for acceptor-FE dissociation to CT
H8 = zeros(num_states);
temp = [2*N_pairs+1, 3*N_pairs+1];
for i=2:2:2*N_pairs
    H8(i, temp) = th;
    temp = temp+N_pairs+1;
    temp= temp(temp<=num_states);
end
H8 = H8+H8';
clear temp
    

H = H1 + H2 + H3 + H4 + H5 + H6 + H7 + H8;
[Psi,D] = eig(H);
% sort eignevalues and eigenvectors
[d,ind] = sort(diag(D),'descend');
D = D(ind,ind);
Psi = Psi(:,ind);

% calculate oscillator strengths
osc_str = abs(sum(Psi)).^2;
for i=1:numel(osc_str)
    osc_str(i) = osc_str(i)*D(i,i);
end

% calculate absorption spectrum
W = 1/pi * gamma./((omega-diag(D)).^2+gamma^2);
abs_spect = osc_str*W;

figure;
subplot(2,2,1)
plot(diag(D), '.')
subplot(2,2,2)
plot(osc_str)
subplot(2,2,3)
plot(1242./omega, sum(W,1))

% % % hold on
% % % plot(omega, sum(W,1));
