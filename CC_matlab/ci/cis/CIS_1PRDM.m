clear all
clc
close all

atoms = {'H','H','H','H'};

atom_valency = [1, 1, 1, 1];

Nelec = sum(atom_valency);

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-7;
ccopts.shift = 0;

lccopts.diis_size = 5;
lccopts.maxit = 150;
lccopts.tol = ccopts.tol;
lccopts.shift = 0.0;

scfopts.diis_size = 5; 
scfopts.maxit = 200; 
scfopts.tol = 1e-8;
scfopts.level_shift = 0.0; 

%
HH_dist = 0.1:0.1:2;
num_pts = length(HH_dist);
COORDS = cell(1,num_pts);

%%
ENERGY_RHF = zeros(1,num_pts);
ENERGY_CC = zeros(1,num_pts);
ENERGY_CRCC = zeros(1,num_pts);

for JJ = 1

        fprintf('\nGEOMETRY %d (OUT OF %d) - \n',JJ,num_pts)
        
        
        COORDS{JJ} = [-HH_dist(JJ)/2, HH_dist(JJ)/2, 0;
                      HH_dist(JJ)/2, HH_dist(JJ)/2, 0;
                      HH_dist(JJ)/2, -HH_dist(JJ)/2, 0;
                     -HH_dist(JJ)/2, -HH_dist(JJ)/2, 0];

        atom_coordinates = COORDS{JJ};

        % get basis set
        basis = get_basis_ccpvdz(atoms, atom_coordinates);

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
            [Escf, C, P, eps_mo, Fock] = rhf(sys_scf,scfopts,P);
        end    

        % record ground state
        ENERGY_RHF(JJ) = Escf;
          
        % ao to mo transformation
        e1int = ao_to_mo(sys_scf.e1int,C);
        e2int = ao_to_mo(sys_scf.e2int,C);
        
        % build up CC system
        nfzc = 0; Nocc_a = Nelec/2; Nocc_b = Nelec/2;
        sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc);

        % UCCSD
        if JJ == 1 % use 0 T vectors as guess
            [cc_t,Ecorr_cc] = uccsd(sys,ccopts);
        else % use previous cluster amplitudes as guess
            [cc_t,Ecorr_cc] = uccsd(sys,ccopts,T_init);
        end
        T_init = cat(1,cc_t.t1a(:),cc_t.t1b(:),cc_t.t2a(:),cc_t.t2b(:),cc_t.t2c(:));

        % UCCSD HBar
        flag_3body = true;
        [HBar_t] = build_ucc_HBar( cc_t, sys, flag_3body);

        % Left UCCSD
        [cc_t,lccsd_resid] = luccsd(cc_t,HBar_t,sys,lccopts);
        
        % CR-CC(2,3)
        [Ecrcc23A,Ecrcc23B,Ecrcc23C,Ecrcc23D] = crcc23_wrap(cc_t,HBar_t,sys);

        % record energies
        ENERGY_CC(JJ) = sys.Escf + Ecorr_cc;
        ENERGY_CRCC(JJ) = Ecrcc23D;

        clear sys cc_t 

end


%% CIS calculation
Nstates = 9;
mult = 1;

% build up CC system
nfzc = 0; Nocc_a = Nelec/2; Nocc_b = Nelec/2;
sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc);

[omega,C1] = cis_spinadapt(Nstates,sys,mult);

% Build CIS 1P RDM
GAMMA_oo = zeros(sys.Nocc_alpha,sys.Nocc_alpha);
GAMMA_vv = zeros(sys.Nvir_alpha,sys.Nvir_alpha);
GAMMA_ov = zeros(sys.Nocc_alpha,sys.Nvir_alpha);
GAMMA_vo = zeros(sys.Nvir_alpha,sys.Nocc_alpha);

for J = 1:Nstates

    c1 = reshape(C1(:,J),[sys.Nvir_alpha,sys.Nocc_alpha]);

    gamma_oo = zeros(sys.Nocc_alpha,sys.Nocc_alpha);
    gamma_vv = zeros(sys.Nvir_alpha,sys.Nvir_alpha);
    gamma_ov = zeros(sys.Nocc_alpha,sys.Nvir_alpha);
    gamma_vo = zeros(sys.Nvir_alpha,sys.Nocc_alpha);

    for i = 1:sys.Nocc_alpha
        for j = 1:sys.Nocc_alpha
            temp = 0.0;
            for a = 1:sys.Nvir_alpha
                temp = temp + 2*c1(a,i)*c1(a,j);
            end
            gamma_oo(i,j) = -temp;
        end
    end


    for a = 1:sys.Nvir_alpha
        for b = 1:sys.Nvir_alpha
            temp = 0.0;
            for i = 1:sys.Nocc_alpha
                temp = temp + 2*c1(a,i)*c1(b,i);
            end
            gamma_vv(a,b) = temp;
        end
    end

    % assuming |CIS> = |HF> + sum_{ai} c1(a,i)*|ia> intermediate
    % normalization. Therefore, 
    % <CIS|N[X^aX_i]|CIS> = sum_{bj} c(b,j)*<HF|N[X^aX_i]|jb> = c(a,i) 
    for a = 1:sys.Nvir_alpha
        for i = 1:sys.Nocc_alpha
            %gamma_ov(i,a) = 2*c1(a,i);
            gamma_ov(i,a) = 0;
        end
    end

    for a = 1:sys.Nvir_alpha
        for i = 1:sys.Nocc_alpha
            %gamma_vo(a,i) = 2*c1(a,i);
            gamma_vo(a,i) = 0;
        end
    end
    

    GAMMA_oo = GAMMA_oo + gamma_oo;
    GAMMA_vv = GAMMA_vv + gamma_vv;
    GAMMA_ov = GAMMA_ov + gamma_ov;
    GAMMA_vo = GAMMA_vo + gamma_vo;

end
% 
% GAMMA_ov = zeros(sys.Nocc_alpha,sys.Nvir_alpha);
% GAMMA_vo = zeros(sys.Nvir_alpha,sys.Nocc_alpha);

GAMMA = [GAMMA_oo, GAMMA_ov;
         GAMMA_vo, GAMMA_vv];

GAMMA = 1/(Nstates+1)*(GAMMA + P);

% Lowdin population matrices
% <psi|P_uv|psi> = sum_{pq} c_p*c_q <p|u>*P_uv*<v|q>
% GAMMA_p = sqrtm(Smat)*GAMMA*sqrtm(Smat);
% P_scf_p = sqrtm(Smat)*P_scf*sqrtm(Smat);
% GAMMA = 1/(Nstates+1)*(GAMMA_p + P_scf_p);

[U,N] = eigenshuffle(GAMMA);


fprintf('MO Label      CIS Occupancy       Occupied\n')
fprintf('---------------------------------------------\n')
for i = 1:sys.Norb
    
    if i <= sys.Nocc_alpha
        fprintf('   %d            %4.8f           *\n',i,N(i))
    else
        fprintf('   %d            %4.8f            \n',i,N(i))
    end

end




