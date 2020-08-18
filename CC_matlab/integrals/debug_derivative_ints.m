clear all
clc
close all

delta = 0.000001;

atoms = {'C','C'};

atom_coordinates = [0, 0, 0;
                    0.5, 0, 0];
                
atom_valency = [6,6];
                
atom_coordinates_dx = atom_coordinates + [0, 0, 0; delta, 0, 0];
                
grad_coords = atom_coordinates(2,:);

basis = get_basis_sto3G(atoms, atom_coordinates);

basis_dx = cell(1,length(basis));
for j = 1:length(basis)
    orb = basis{j};
    if orb.origin == grad_coords
        orb.origin = orb.origin + [delta, 0, 0];
    end
    basis_dx{j} = orb;
end

%% Testing derivative integrals

Smat = get_Smat(basis); Smat_dx = get_Smat(basis_dx);
Zmat = get_Zmat(basis,atom_coordinates,atom_valency); Zmat_dx = get_Zmat(basis_dx,atom_coordinates_dx,atom_valency);
VVmat = get_ERImat(basis,0.0); VVmat_dx = get_ERImat(basis_dx,0.0);

dSmat_numerical = 1/delta.*(Smat_dx-Smat);
dZmat_numerical = 1/delta.*(Zmat_dx-Zmat);
dVVmat_numerical = 1/delta.*(VVmat_dx-VVmat);

dSmat = get_dSmat(basis,grad_coords);
get_error(squeeze(dSmat{1}(:,:,1)),dSmat_numerical)
dSj = squeeze(dSmat{1}(:,:,1));

dZmat = get_dZmat(basis,atom_coordinates,atom_valency,grad_coords);
get_error(squeeze(dZmat{1}(:,:,1)),dZmat_numerical)
dZj = squeeze(dZmat{1}(:,:,1));

dVVmat = get_dERImat(basis,grad_coords);
get_error(squeeze(dVVmat{1}(:,:,:,:,1)),dVVmat_numerical)
dVVj = squeeze(dVVmat{1}(:,:,:,:,1));






