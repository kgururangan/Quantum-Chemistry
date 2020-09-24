clear all
clc
close all

% [TODO Wishlist]

% Coupled-cluster
% (1) Add t1b and t2c updates in UCCSD
% (2) Spin-integrated HBar diagonal (S, D, and T diagonals for H1, H2, H3)
% (3) Spin-integrated CR-CC(2,3) and CR-EOMCC(2,3)
% (4) Spin-integrated CCSDt(II) and CCSDt(III)
% (5) Spin-integrated CCSDt(II) HBar
% (6) Spin-integrated left-CCSDt(II)
% (7) Spin-integrated CC(t,3)(II) using loop over full triples space
%     and just accumulating correction corresponding to amplitudes in
%     Q-space
% (8) Spin-integrated EOMCCSDt(II)

% SCF
% (1) Implement CPHF equation solver - analytic gradients and hessians
% (2) UHF and GHF solvers for spin-frustrated systems
% (3) Relativistic Dirac-Hartree-FOck solver
% (4) MCSCF and CASSCF solvers

% Integrals and SCF
% (1) Symmetry-adapt AO basis function to D2h and subgroups
% (2) Replace MMD integral scheme with Obara-Saika
% (3) Implement density fitting for ERI
% (4) Derivative integrals using MMD or Obara-Saika


%%
flag_check = false;

if flag_check
    load h2o-pvdz-ao-ints.mat
end
%
Nelec = 12;

atoms = {'C','C'};

atom_valency = [6, 6];

geometry = '2.5Re';

switch geometry
    case 'Re'
        atom_coordinates = [0, 1.515, -1.058;
                            0, -1.515, -1.058;
                            0, 0, -0.0090];
    case '1.5Re'
        atom_coordinates = [0, 2.272, -1.588;
                            0, -2.272, -1.588;
                            0, 0, -0.0135];
    case '2Re'
        atom_coordinates = [0, 3.030, -2.117;
                            0, -3.030, -2.117;
                            0, 0, -0.0180];
    case '2.5Re'
        atom_coordinates = [0, 3.788, -2.647;
                            0, -3.788, -2.647;
                            0, 0, -0.0225];
    case '3Re'
        atom_coordinates = [0, 4.545, -3.176;
                            0, -4.545, -3.176;
                            0, 0, -0.0270];
    case '8Re'
        atom_coordinates = [0, 12.122, -8.471;
                            0, -12.122, -8.471;
                            0, 0, -0.072];
    otherwise
        disp('Please select valid H2O geometry from Bartlett, Balkova paper!')
end

R0 = 1.3983972315484032;
HOH_angle = deg2rad(114.5);
theta = 0.5*(pi - HOH_angle);
atom_coordinates = [ R0*cos(theta), R0*sin(theta), 0;
                    -R0*cos(theta), R0*sin(theta), 0;
                          0       ,      0       , 0 ];  

atom_coordinates = [0, 0, -0.739;
                    0, 0, 0.739];

basis = get_basis_ccpvdz(atoms, atom_coordinates);

% gb1.shell = [0, 0, 0];
% gb1.coef = [8.048803E-06, ...
%           6.258306E-05, ...
%           3.290239E-04, ...
%            1.387355E-03, ...
%            5.023256E-03, ...
%           1.610140E-02, ...
%           4.590034E-02, ...
%           1.136154E-01, ...
%            2.283869E-01, ...
%           3.221159E-01, ...
%            2.383661E-01 , ...
%           7.404667E-02, ...
%           9.214197E-02, ...
%            9.339790E-02, ...
%           1.573965E-02, ...
%          -4.186682E-04, ...
%            5.376318E-05, ...
%           -3.816654E-05 ...
%          4.319603E-05, ...
%          -3.401019E-06];
%          
% gb1.exps = [4.316265E+06,  ...       
%          6.463424E+05    ,   ...
%          1.470897E+05     ,   ...
%          4.166152E+04     ,    ...
%          1.359077E+04      ,   ...
%          4.905750E+03      ,...
%          1.912746E+03       ,  ... 
%          7.926043E+02      ,  ...
%          3.448065E+02      , ...
%         1.558999E+02        ,...
%         7.223091E+01    ,...
%         3.272506E+01     , ...  
%         1.566762E+01      ,  ...
%         7.503483E+00       ,  ... 
%         3.312223E+00        ,  ...
%         1.558471E+00        ,...
%         6.839140E-01      ,...
%         1.467570E-01       , ... 
%         7.058300E-02       ,...
%         3.144900E-02];      
%  gb1.origin = atom_coordinates(1,:);
%  gb2 = gb1; gb2.origin = atom_coordinates(2,:);
% 
%     
% basis = {gb1, gb2};
% 
Norb = length(basis);
% 
% mol.Nelec = Nelec;
% mol.atoms = atoms;
% mol.atom_coordinates = atom_coordinates;
% mol.atom_valency = atom_valency;
% mol.basis = basis;

Vnuc = calc_nuclear_nuclear(atom_coordinates,atom_valency);

sys.Vnuc = Vnuc;


%%
tic
Smat = get_Smat(basis);
toc

if flag_check
    errS = Smat - Smat_py;
    sum(abs(errS(:)))
end
%%

tic
Zmat = get_Zmat(basis,atom_coordinates,atom_valency);
toc

if flag_check
    errZ = Zmat - Zmat_py;
    sum(abs(errZ(:)))
end
%%

[VVmat] = get_ERImat(basis,0.0);

if flag_check
    errVV = VVmat - permute(VVmat_py,[1,3,2,4]);
    sum(abs(errVV(:)))
end

%%

sys.e1int = Zmat; sys.overlap = Smat; sys.e2int = permute(VVmat,[1,3,2,4]); 
sys.Vnuc = Vnuc; sys.Nelec = 10;

opts.diis_size = 5; opts.maxit = 200; opts.tol = 1e-8;
opts.level_shift = 0.0;

[Escf, C, P, eps_mo, Fock] = rhf(sys,opts);

%% Test H2 cc-pvdz integrals
clear all
clc
close all

atoms = {'H', 'H'};
atom_valency = [1, 1];
Nelec = 2;

atom_coordinates = [0, 0, -0.6991986157742016;
                    0, 0,  0.6991986157742016]; 
                
basis = get_basis_ccpvdz(atoms, atom_coordinates);

Norb = length(basis);

Smat = get_Smat(basis);
Zmat = get_Zmat(basis,atom_coordinates,atom_valency);
[VVmat] = get_ERImat(basis,0.0);

% for i = 1:Norb
%     for j = i:Norb
%         for k = j:Norb
%             for l = k:Norb
%                 fprintf('V(%d,%d,%d,%d) = %4.10f\n',i,j,k,l,VVmat(i,j,k,l))
%             end
%         end
%     end
% end
%%

err = zeros(Norb,Norb,Norb,Norb);
idx = [];
ct = 1;
for a = 1:Norb
    for b = 1:Norb
        for c = 1:Norb
            for d = 1:Norb
                err(a,b,c,d) = VVmat_py(a,c,b,d) - VVmat(a,b,c,d);
                if abs(err(a,b,c,d)) > 1e-6
                    idx(ct,:) = [a,b,c,d];
                    ct = ct + 1;
                end
            end
        end
    end
end
%%
errf = @(m) err(idx(m,1),idx(m,2),idx(m,3),idx(m,4))

%%

n = 0:8;
T = linspace(0,25,500);
f1 = zeros(length(n),length(T));
for j = 1:length(n)
    for i = 1:length(T)
        f1(j,i) = boys(n(j),T(i));
    end
end


