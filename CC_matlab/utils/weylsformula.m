%%
clc
clear all
close all

zH = 1;
zC = 6;
zO = 8;

h_2z = 5;
h_3z = 14;
h_4z = 30;

o_2z = 14;
o_3z = 30;
o_4z = 55;

c_2z = 14;
c_3z = 30;
c_4z = 55;

S = 0;

%% C2H4

Ntot = 2*zC + 2*zH;
Nfz = 4;
N = Ntot - Nfz;
fprintf('===================++C2H2 (N=%d, S=%d)++===================\n',N,S)
fprintf('Basis          N(CSF)              K(Orbitals)\n')
% pvdz
K = 2*h_2z + 2*c_2z - Nfz/2;
fprintf('cc-pVDZ        %4.6e        %d\n',weyl(K,N,S),K)
% pvtz
K = 2*h_3z + 2*c_3z - Nfz/2;
fprintf('cc-pVTZ        %4.6e        %d\n',weyl(K,N,S),K)
% pvqz
K = 2*h_4z + 2*c_4z - Nfz/2;
fprintf('cc-pVQZ        %4.6e        %d\n',weyl(K,N,S),K)

%% C2H4

Ntot = 2*zC + 4*zH;
Nfz = 4;
N = Ntot - Nfz;
fprintf('===================++C2H4 (N=%d, S=%d)++===================\n',N,S)
fprintf('Basis          N(CSF)              K(Orbitals)\n')
% pvdz
K = 2*c_2z + 4*h_2z - Nfz/2;
fprintf('cc-pVDZ        %4.6e        %d\n',weyl(K,N,S),K)
% pvtz
K = 2*c_3z + 4*h_3z - Nfz/2;
fprintf('cc-pVTZ        %4.6e        %d\n',weyl(K,N,S),K)
% pvqz
K = 2*c_4z + 4*h_4z - Nfz/2;
fprintf('cc-pVQZ        %4.6e        %d\n',weyl(K,N,S),K)

%% C6H6

Ntot = 6*zC + 6*zH;
Nfz = 12;
N = Ntot - Nfz;
fprintf('===================++C6H6 (N=%d, S=%d)++===================\n',N,S)
fprintf('Basis          N(CSF)              K(Orbitals)\n')
% pvdz
K = 6*h_2z + 6*c_2z - Nfz/2;
fprintf('cc-pVDZ        %4.6e        %d\n',weyl(K,N,S),K)
% pvtz
K = 6*h_3z + 6*c_3z - Nfz/2;
fprintf('cc-pVTZ        %4.6e        %d\n',weyl(K,N,S),K)
% pvqz
K = 6*h_4z + 6*c_4z - Nfz/2;
fprintf('cc-pVQZ        %4.6e        %d\n',weyl(K,N,S),K)

%% H4
Ntot = 4*zH;
Nfz = 0;
N = Ntot - Nfz;
K = 4*47 - Nfz/2;
S = 0;
fprintf('aug-cc-pVQZ        %4.6e        %d\n',weyl(K,N,S),K)

%% C4H4 / cc-pVDZ

Ntot = 4*zH + 4*zC;
Nfz = 8;
N = Ntot - Nfz;
K = 76 - Nfz/2;
S = 0;
fprintf('cc-pVDZ        %4.6e        %d\n',weyl(K,N,S),K)

%% F2 / cc-pVDZ

Ntot = 18;
Nfz = 4;
N = Ntot - Nfz;
K = 30 - Nfz/2;
S = 0;
fprintf('cc-pVDZ        %4.6e        %d\n',weyl(K,N,S),K)

%% F2 / cc-pVTZ

Ntot = 18;
Nfz = 4;
N = Ntot - Nfz;
K = 60 - Nfz/2;
S = 0;
fprintf('cc-pVDZ        %4.6e        %d\n',weyl(K,N,S),K)

%%
function [csf] = weyl(K,N,S)

%     K = 264;
%     N = 30;
%     S = 0;

    % K+1 choose N/2-S

    csf = (2*S+1)/(K+1);

    c1 = 1;
    for j = K+1-(N/2-S)+1 : K+1
        c1 = c1 * j;
    end

    c1;

    csf = csf * c1 / factorial(N/2-S);

    % K+1 choose N/2+S+1
    c2 = 1;
    for j = K+1-(N/2+S+1)+1 : K+1
        c2 = c2 * j;
    end

    c2;

    csf = csf * c2/factorial(N/2+S+1);

end


