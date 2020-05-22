clear all
clc
close all

addpath(genpath('C:\Users\Karthik\Dropbox\Hartree Fock\hartree_fock\CT_Hamiltonian\Dimer_FCI'))
load monomer_ao_mo.mat

Norb = size(xi_mat,1);

atom_coords = unique(origin_mat,'rows','stable');

Natom = size(atom_coords,1);

delta = 6;
R_monomer_1 = [0,0,0];
R_monomer_2 = [0,0,-delta/2];

coords_1 = zeros(Natom,3); coords_2 = zeros(Natom,3);

for i = 1:Natom
    coords_1(i,:) = atom_coords(i,:)-R_monomer_1;
    coords_2(i,:) = atom_coords(i,:)-R_monomer_2;
end

orbs_fcn_1 = cell(1,Norb);
orbs_fcn_2 = cell(1,Norb);
for i = 1:Norb
    orbs_fcn_1{i} = ao_fcn(shell_mat(i,:),xi_mat(i,:),coeff_mat(i,:),origin_mat(i,:)-R_monomer_1);
    orbs_fcn_2{i} = ao_fcn(shell_mat(i,:),xi_mat(i,:),coeff_mat(i,:),origin_mat(i,:)-R_monomer_2);
end
%
    
% bohr
dist_x = max(origin_mat(:,1)) - min(origin_mat(:,1));
dist_y = max(origin_mat(:,2)) - min(origin_mat(:,2));
dist_z = max(origin_mat(:,3)) - min(origin_mat(:,3));

fact = 5;

Lx = fact*dist_x;
Ly = fact*dist_y;
%Lz = fact*dist_z;
Lz = 40;

h = 0.8;

dx = h; dy = h; dz = h;

x = -Lx/2:dx:Lx/2; Nx = length(x);
y = -Ly/2:dy:Ly/2; Ny = length(y);
z = -Lz/2:dz:Lz/2; Nz = length(z);

[X,Y,Z] = meshgrid(x,y,z); % meshgrid outpute length(y) x length(x) x length(z)
X = permute(X,[2,1,3]); Y = permute(Y,[2,1,3]); Z = permute(Z,[2,1,3]);

LUMO1 = zeros(Nx,Ny,Nz); LUMO2 = zeros(Nx,Ny,Nz);
HOMO1 = zeros(Nx,Ny,Nz); HOMO2 = zeros(Nx,Ny,Nz);

for J = 1:Norb
    J
    
    LUMO1 = LUMO1 + lumo_vec(J)*orbs_fcn_1{J}(X,Y,Z);
    HOMO1 = HOMO1 + homo_vec(J)*orbs_fcn_1{J}(X,Y,Z);
    
    LUMO2 = LUMO2 + lumo_vec(J)*orbs_fcn_2{J}(X,Y,Z);
    HOMO2 = HOMO2 + homo_vec(J)*orbs_fcn_2{J}(X,Y,Z);
end
%
ORBS = {HOMO1, HOMO2, LUMO1, LUMO2};

%%

% MO vectors not orthonormal among themselves! they are orthonormal wrt to S (generalized
% eigenvalue problem)... also contracted gaussians do not necesarilly
% integrate to 1

% The overlap of MO's should be identity though...
% are we loading things from the gamess log file incorrectly?


Smat = zeros(4);
for i = 1:4
    for j = 1:4
        Smat(i,j) = sum(sum(sum(ORBS{i}.*ORBS{j},1)*dx,2)*dy,3)*dz;
    end
end
Smat
Xorth = orthomat(Smat,'symmetric');

%% Plotting HOMO and LUMO in 3D

figure(1)

subplot(121)
scatter3(coords_1(:,1),coords_1(:,2),coords_1(:,3),40,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0]); 
hold on 
% scatter3(coords_2(:,1),coords_2(:,2),coords_2(:,3),40,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1]); 
%
C = parula(12036);
p = patch(isosurface(X,Y,Z,HOMO1,'noshare'),'FaceVertexCData',C);
    length(p.Vertices) % check this to set size of C
    p.FaceColor = 'none';
    p.EdgeColor = 'interp';
    view(50,40)
    
% C = parula(10572);
% p = patch(isosurface(X,Y,Z,HOMO2,'noshare'),'FaceVertexCData',C);
%     length(p.Vertices) % check this to set size of C
%     p.FaceColor = 'none';
%     p.EdgeColor = 'interp';
%     view(50,40)
%       
hold off
xlabel('y / bohr'); ylabel('x / bohr'); zlabel('z / bohr');
set(gca,'FontSize',14,'Linewidth',2,'Box','off')
axis square
%axis([x(1),x(end),y(1),y(end),z(1),z(end)])
title('HOMO')

subplot(122)
scatter3(coords_1(:,1),coords_1(:,2),coords_1(:,3),40,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0]); 
hold on 
% scatter3(coords_2(:,1),coords_2(:,2),coords_2(:,3),40,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1]); 

C = parula(7524);
p = patch(isosurface(X,Y,Z,LUMO1,'noshare'),'FaceVertexCData',C);
    length(p.Vertices) % check this to set size of C
    p.FaceColor = 'none';
    p.EdgeColor = 'interp';
    view(50,40)
    
% C = parula(24036);
% p = patch(isosurface(X,Y,Z,LUMO2,'noshare'),'FaceVertexCData',C);
%     length(p.Vertices) % check this to set size of C
%     p.FaceColor = 'none';
%     p.EdgeColor = 'interp';
%     view(50,40) 
    
hold off
xlabel('y / bohr'); ylabel('x / bohr'); zlabel('z / bohr');
set(gca,'FontSize',14,'Linewidth',2,'Box','off')
axis square
title('LUMO')


%% Compute ERI

lambda = 0;
[VVmat] = eri_fft({HOMO1,HOMO2,LUMO1,LUMO2},[Lx,Ly,Lz],[dx,dy,dz],lambda)

%%


function [orb] = ao_fcn(shell,xi,coeff,origin)

        L = shell(1);
        M = shell(2);
        N = shell(3);
        
        % AOs need primitive normalization factor        
        orb = @(x,y,z) (x-origin(1)).^L.*(y-origin(2)).^M.*(z-origin(3)).^N.*(...
                        coeff(1)*primitive_norm(L,M,N,xi(1))*exp(-xi(1)*((x-origin(1)).^2+(y-origin(2)).^2+(z-origin(3)).^2))+...
                        coeff(2)*primitive_norm(L,M,N,xi(2))*exp(-xi(2)*((x-origin(1)).^2+(y-origin(2)).^2+(z-origin(3)).^2))+...
                        coeff(3)*primitive_norm(L,M,N,xi(3))*exp(-xi(3)*((x-origin(1)).^2+(y-origin(2)).^2+(z-origin(3)).^2)) );

        %         
        % orb = @(x,y,z) (x-origin(1)).^L.*(y-origin(2)).^M.*(z-origin(3)).^N.*(...
        %                 coeff(1)*exp(-xi(1)*((x-origin(1)).^2+(y-origin(2)).^2+(z-origin(3)).^2))+...
        %                 coeff(2)*exp(-xi(2)*((x-origin(1)).^2+(y-origin(2)).^2+(z-origin(3)).^2))+...
        %                 coeff(3)*exp(-xi(3)*((x-origin(1)).^2+(y-origin(2)).^2+(z-origin(3)).^2)) );

end
        