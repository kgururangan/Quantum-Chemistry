clear all
clc
close all
% 
% fid = fopen('monomer_geometry.inp');
% 
% atoms = {};
% atom_valency = [];
% atom_coordinates = [];
% ct = 1;
% 
% tline = fgetl(fid);
% while ischar(tline)
%     L = split(tline);
%     if length(L) > 1
%         atoms{ct} = L{1};
%         atom_valency(ct) = str2double(L{2});
%         atom_coordinates(ct,:) = [str2double(L{3}), str2double(L{4}), str2double(L{5})];
%         ct = ct + 1;
%     end
%     %disp(tline)
%     tline = fgetl(fid);
% end

%

load('/Users/karthik/Desktop/pent_monomer_SF/pent_monomer_geometry.mat')

%



scfopts.diis_size = 3; 
scfopts.maxit = 5; 
scfopts.tol = 1e-8;
scfopts.level_shift = 0.0; 

Nelec = sum(atom_valency);
Nocc = Nelec/2;

% get basis set
basis = get_basis_sto3G(atoms, atom_coordinates);
Norb = length(basis);

% calculate AO integrals
Vnuc = calc_nuclear_nuclear(atom_coordinates,atom_valency);
Smat = get_Smat(basis); 
Zmat = get_Zmat(basis,atom_coordinates,atom_valency); 
[VVmat_chemist] = get_ERImat(basis,0.0);

% SCF solver
sys_scf.e1int = Zmat; sys_scf.overlap = Smat; 
sys_scf.e2int = permute(VVmat_chemist,[1,3,2,4]);
sys_scf.Vnuc = Vnuc; sys_scf.Nelec = Nelec; 
if JJ == 1
    [Escf, C, P, eps_mo, Fock] = rhf(sys_scf,scfopts);
else
    [Escf, C, P, eps_mo, Fock] = rhf(sys_scf,scfopts);
end
