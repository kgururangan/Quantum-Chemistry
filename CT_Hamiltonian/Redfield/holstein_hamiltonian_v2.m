function [Htot,FCmat] = holstein_hamiltonian_v2(n_quanta,omega_g,omega_ge,numstates,omega_e,Dpv,Jcoul,G,zp)

n_modes = length(omega_g);
N = n_quanta^n_modes;
Nproblem = numstates*N;

% replaced Bv with omega_e
Bv = sqrt(omega_e./omega_g);
Bv = [ones(1,n_modes);Bv];


GmnE = cell(numstates,numstates);
H = cell(numstates,numstates);

% Step 1: define the SHO vibrational basis operators
Cre = zeros(n_quanta,n_quanta);
Ann = zeros(n_quanta,n_quanta);

for i = 0:n_quanta-1
    for j = 0:n_quanta-1
        
        if j == i + 1
            Ann(i+1,j+1) = sqrt(j);
        end
        
        if i == j + 1
            Cre(i+1,j+1) = sqrt(i);
        end
    end
end
Num = Cre*Ann;

% Step 2: Full vibronic Hamiltonian matrix
for k = 1:numstates
    Htemp = zeros(N,N);
    for i = 1:n_modes
        if zp ~= 0
            H0temp = omega_g(i)*(Num + 0.5*eye(n_quanta));
        else
            H0temp = omega_e(k,i)*Num + Dpv(k,i)*omega_e(k,i)*(Cre+Ann);% + Dpv(k,i)^2*Bv(k,i)*omega_g(k,i);
        end
        A = 1;
        for j = 1:n_modes
            if j ~= i
                A = kron(eye(n_quanta),A);
            else
                A = kron(H0temp,A);
            end
        end
        Htemp = Htemp + A;
    end
    
    GmnE{1,k} = myFranckCondon(Bv(k,:), Dpv(k,:), n_modes, n_quanta);
    GmnE{k,1} = GmnE{1,k}';
    
    GmnE{k,k} = zeros(N);
    H{k,k} = Htemp + omega_ge(k)*eye(N);
end


for k = 1:numstates
    for kp = 1:numstates
        if k ~= kp
            if k > 1 && kp > 1
                GmnE{k,kp} = transpose(GmnE{1,k})*GmnE{1,kp};
            end
            H{k,kp} = Jcoul(k,kp)*eye(N);
        end
    end
end

% Anharmonic correction to ground state
Hanh_g = zeros(N,N); % ground-state harmonic Hamiltonian

% Hanh_g = kron(kron(Ann + Cre,Ann + Cre),Ann + Cre);
for i = 1:n_modes
    for j = 1:n_modes
        for k = 1:n_modes
            vec = [i,j,k];
            ALL = 1;
            for l = 1:n_modes
               [Il,Jl] = find(vec == l);
               nl = length(Jl);
               ALL = kron((Ann + Cre)^nl,ALL);
            end
            Hanh_g = Hanh_g + G(i,j,k)*ALL;
        end
    end
end

Hanh_g = (1/6)*Hanh_g;
H{1,1} = H{1,1} + Hanh_g;
%Htot(1:N,1:N) = Htot(1:N,1:N) + Hanh_g;

Htot = cell2mat(H);
FCmat = cell2mat(GmnE);

end