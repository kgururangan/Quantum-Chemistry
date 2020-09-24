clear all
clc
close all

atoms = {'H','H'};

atom_valency = [1, 1];

Natoms = length(atoms);
Nelec = sum(atom_valency);
Nocc = Nelec/2;

atom_coordinates = [0, 0, 0;
                   0, 0, 0.1];
               
x0 = atom_coordinates(:);

v1 = zeros(length(x0),1); v1(1:Natoms:end) = 1;
v2 = zeros(length(x0),1); v2(2:Natoms:end) = 1;
v3 = zeros(length(x0),1); v3(3:Natoms:end) = 1;

posv = {v1, v2, v3};

Hessian = eye(size(x0,1));
L = chol(Hessian,'lower');

scf_out = scf_fcn(x0, atoms, atom_valency);

Escf = scf_out{1};
MO_vec = scf_out{2};
eps_mo = scf_out{3};

[dE_rhf0] = gradscf_fcn(x0,atoms,atom_valency,MO_vec,eps_mo);

%%
%
it = 0; maxit = 1000;
while it <= maxit && norm(dE_rhf0) > 1e-6

        fprintf('\nGEOMETRY %d - \n',it)
        
        search_dir = -L'\(L\dE_rhf0);
        
        fprintf('\nBeginning Golden Section line search algorithm...\n')
        
        fun = @(alpha) scf_fcn(x0+alpha*search_dir, atoms, atom_valency);
        
        a = 0; b = 5;                   % search interval range        
        tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
        k=0;                            % number of iterations
        iter = 25;                      % max no. iterations
        epsilon = 1e-6;                 % interval tolerance

        x1=a+(1-tau)*(b-a);             % computing x values
        x2=a+tau*(b-a);
        
        % computing values in x points
        f_x1 = fun(x1);       
        Escf1 = f_x1{1}; C1 = f_x1{2};
        
        f_x2 = fun(x2);
        Escf2 = f_x2{1}; C2 = f_x2{2};

        while ((abs(b-a)>epsilon) && (k<iter))
            
            fprintf('GOLDEN SECTION IT - %d    INTERVAL = %4.4f\n',k,abs(b-a))
            
            k=k+1;
            
            if(Escf1 < Escf2)
                
                b=x2;
                x2=x1;
                x1=a+(1-tau)*(b-a);

                f_x1 = fun(x1);
                Escf1 = f_x1{1}; C1 = f_x1{2}; eps_mo1 = f_x1{3};
                
                f_x2 = fun(x2);
                Escf2 = f_x2{1}; C2 = f_x2{2}; eps_mo2 = f_x2{3};

            else
                
                a=x1;
                x1=x2;
                x2=a+tau*(b-a);

                f_x1 = fun(x1);
                Escf1 = f_x1{1}; C1 = f_x1{2}; eps_mo1 = f_x1{3};
                
                f_x2 = fun(x2);
                Escf2 = f_x2{1}; C2 = f_x2{2}; eps_mo2 = f_x2{3};
                
            end

            k=k+1;
            
        end
        
        fprintf('\nGolden section line search converged...\n')

        % choose minimum point
        if(Escf1 < Escf2)
            alpha = x1;
            f_x1 = fun(x1);
            Escf = f_x1{1}; MO_vec = f_x1{2}; eps_mo = f_x1{3};
        else
            alpha = x2;
            fval = f_x2;
            Escf = f_x2{1}; MO_vec = f_x2{2}; eps_mo = f_x2{3};
        end

        x = x0 + alpha*search_dir;
        
        dE_rhf = gradscf_fcn(x,atoms,atom_valency,MO_vec,eps_mo);
        
        grad_conv = dE_rhf - dE_rhf0;
        dx = x - x0;

        D = (grad_conv*grad_conv')/(grad_conv'*dx);
        E = (dE_rhf0*dE_rhf0')/(dE_rhf0'*search_dir);
        Hessian = Hessian + D + E;
        L = chol(Hessian,'lower');

        dE_rhf0 = dE_rhf;
        x0 = x;

        fprintf('Iter - %d:      FVAL = %4.8f      OPTIM = %4.8f\n',it, Escf, norm(dE_rhf))

        it = it + 1;

end
    
xsol = x;

if norm(dE_rhf) <= 1e-6
    fprintf('BFGS successfuly converged in %d iterations\n',it-1)
else
    fprintf('BFGS not converged\n')
end
        




%%
function [xout] = scf_fcn(x, atoms, atom_valency, scfopts)

    if nargin < 4
        % scf options
        scfopts.diis_size = 3; 
        scfopts.maxit = 200; 
        scfopts.tol = 1e-8;
        scfopts.level_shift = 0.0; 
    end
    
    Nelec = sum(atom_valency); Nocc = Nelec/2;
    
    atom_coordinates = reshape(x,[length(atoms),3]);

    % get basis set
    basis = get_basis_sto3G(atoms, atom_coordinates);

    % calculate AO integrals
    Vnuc = calc_nuclear_nuclear(atom_coordinates,atom_valency);
    Smat = get_Smat(basis); 
    Zmat = get_Zmat(basis,atom_coordinates,atom_valency); 
    [VVmat_chemist] = get_ERImat(basis,0.0);

    % SCF solver
    sys_scf.e1int = Zmat; sys_scf.overlap = Smat; 
    sys_scf.e2int = permute(VVmat_chemist,[1,3,2,4]);
    sys_scf.Vnuc = Vnuc; sys_scf.Nelec = Nelec; 

    [Escf, C, P, eps_mo, Fock] = rhf(sys_scf,scfopts);
    
    xout = {Escf, C, eps_mo};
    
end

function [dE_RHF] = gradscf_fcn(x,atoms,atom_valency,C,eps_mo)

    Nelec = sum(atom_valency); Nocc = Nelec/2;
    
    atom_coordinates = reshape(x,[length(atoms),3]);
    
    basis = get_basis_sto3G(atoms, atom_coordinates);

    % calculate derivative AO integrals
    dSmat = get_dSmat(basis,atom_coordinates);
    dZmat = get_dZmat(basis,atom_coordinates,atom_valency,atom_coordinates);
    dVVmat = get_dERImat(basis,atom_coordinates);
    dVnn = grad_nuclear_nuclear(atom_coordinates,atom_valency,atom_coordinates);
    
    % density matrix
    P = 2*einsum_kg(C(:,1:Nocc),C(:,1:Nocc),'pi,qi->pq');

    % energy-weighted density matrix
    W = 2*einsum_kg(C(:,1:Nocc)*diag(eps_mo(1:Nocc)),C(:,1:Nocc),'pi,qi->pq');

    % calcualte analytical gradient
    dE_RHF = zeros(size(atom_coordinates,1),3);
    for M = 1:size(atom_coordinates,1)
        dZ = dZmat{M}; dS = dSmat{M}; dERI = dVVmat{M}; dVnuc = dVnn{M};
        for j = 1:3
            dZj = squeeze(dZ(:,:,j)); dSj = squeeze(dS(:,:,j)); dERIj = squeeze(dERI(:,:,:,:,j));
            dERIj_asym = dERIj - 0.5*permute(dERIj,[1,4,3,2]);
            dE_RHF(M,j) = einsum_kg(P,dZj,'uv,uv->')...
                             +0.5*einsum_kg(einsum_kg(P,dERIj_asym,'ls,uvls->uv'),P,'uv,uv->')...
                             -einsum_kg(W,dSj,'uv,uv->')...
                             +dVnuc(j);      
        end
    end      

    dE_RHF = dE_RHF(:);
    
end