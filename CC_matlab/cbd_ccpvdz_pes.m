clear all
clc
close all

Nelec = 28;

atoms = {'H','H', 'H', 'H', 'C', 'C', 'C', 'C'};

atom_valency = [1, 1, 1, 1, 6, 6, 6, 6];

nstates = 1; % including ground state

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-10;
ccopts.shift = 0;

lccopts.diis_size = 5;
lccopts.maxit = 150;
lccopts.tol = ccopts.tol;
lccopts.shift = 0.0;

eomopts.nroot = nstates-1;
eomopts.maxit = 300;
eomopts.tol = 1e-7;
eomopts.nvec_per_root = 1;
eomopts.max_nvec_per_root = 5;
eomopts.flag_verbose = 1;
eomopts.mult = 1;
eomopts.thresh_vec = 1e-3;
eomopts.solver = 2;

leomccopts.maxit = 150;
leomccopts.diis_size = 5;
leomccopts.tol = 1e-6; 
leomccopts.shift = 0;
leomccopts.solver = 1;
leomccopts.nroot = nstates-1;

scfopts.diis_size = 5; 
scfopts.maxit = 200; 
scfopts.tol = 1e-8;
scfopts.level_shift = 0.0; 

%

theta = linspace(20,45,20);
CC_dist = linspace(1,4,20);
num_pts = length(CC_dist); 
COORDS = cell(1,num_pts);
for k = 1:num_pts
    at_H1 = [0, 1.55801763, 0];
    at_H2 = [0, 1.55801763, 0];
    at_O = [0, 0, 0];
    COORDS{k} = [at_H1; at_H2; at_O];
end

%%
ENERGY_RHF = zeros(1,num_pts);
ENERGY_CC = zeros(nstates,num_pts);
ENERGY_CRCC = zeros(nstates,num_pts);


int_lab = {'1.0','1.1','1.2','1.3','1.4','1.5','1.6','1.7','1.8','1.9','2.05' ...
           '2.16','2.26','2.37','2.47','2.58','2.68','2.79','2.89','3.0'};

for JJ = 1:num_pts

        fprintf('\nGEOMETRY %d (OUT OF %d) - \n',JJ,num_pts)

%         atom_coordinates = COORDS{JJ};
% 
%         % get basis set
%         basis = get_basis_ccpvdz(atoms, atom_coordinates);
% 
%         % calculate AO integrals
%         Vnuc = calc_nuclear_nuclear(atom_coordinates,atom_valency);
%         Smat = get_Smat(basis); 
%         Zmat = get_Zmat(basis,atom_coordinates,atom_valency); 
%         [VVmat_chemist] = get_ERImat(basis,0.0);
%       
%         % SCF solver
%         sys_scf.e1int = Zmat; sys_scf.overlap = Smat; 
%         sys_scf.e2int = permute(VVmat_chemist,[1,3,2,4]);
%         sys_scf.Vnuc = Vnuc; sys_scf.Nelec = 10; 
%   
%         if JJ == 1
%             [Escf, C, P, eps_mo, Fock] = rhf(sys_scf,scfopts);
%         else
%             [Escf, C, P, eps_mo, Fock] = rhf(sys_scf,scfopts,P);
%         end
%
%         % record ground state
%         ENERGY_RHF(JJ) = sys.Escf;
%           
%         % ao to mo transformation
%         e1int = ao_to_mo(sys_scf.e1int,C);
%         e2int = ao_to_mo(sys_scf.e2int,C);

        % load MO integrals from GAMESS
        workpath = '/Users/harellab/Desktop/CC_matlab_tests/h2o-pvdz-gms/PES/';
        workpath = strcat(strcat(strcat(workpath,'h2o-'),int_lab{JJ}),'/');
        [e1int, e2int, Vnuc, Norb] = load_integrals(workpath);

        % build up CC system
        nfzc = 0; Nocc_a = Nelec/2; Nocc_b = Nelec/2;
        sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc);

        % UCCSD
        if JJ == 1 % use 0 T vectors as guess
            [cc_t,Ecorr_ucc] = uccsd(sys,ccopts);
        else % use previous cluster amplitudes as guess
            [cc_t,Ecorr_cc] = uccsd(sys,ccopts,T_init);
        end
        T_init = cat(1,cc_t.t1a(:),cc_t.t1b(:),cc_t.t2a(:),cc_t.t2b(:),cc_t.t2c(:));

        % UCCSD HBar
        flag_3body = true;
        [HBar_t] = build_ucc_HBar( cc_t, sys, flag_3body);

        % Left UCCSD
        [cc_t,lccsd_resid] = luccsd(cc_t,HBar_t,sys,lccopts);

        % EOM-UCCSD
        %if JJ == 1 % use CIS guess
            eomopts.init_guess = 'cis';
            [Rvec, omega, eom_residual, cc_t] = eomuccsd(HBar_t,cc_t,sys,eomopts);
        %else % use previous R vectors as guess
        %    eomopts.init_guess = 'custom';
        %    [Rvec, omega, eom_residual, cc_t] = eomuccsd(HBar_t,cc_t,sys,eomopts,Rvec);
        %end

        % Left EOM-UCCSD
        [Lvec,eom_lcc_resid,cc_t] = lefteomuccsd(omega,Rvec,HBar_t,cc_t,sys,leomccopts);

        % CR-CC(2,3)
        [Ecrcc23A,Ecrcc23B,Ecrcc23C,Ecrcc23D] = crcc23_wrap(cc_t,HBar_t,sys,omega);

        % record energies
        for n = 1:nstates
            if n == 1
                ENERGY_CC(n,JJ) = sys.Escf + Ecorr_ucc;
            else
                ENERGY_CC(n,JJ) = sys.Escf + Ecorr_ucc + omega(n-1);
            end
            ENERGY_CRCC(n,JJ) = Ecrcc23D(n);
        end

        clear sys cc_t 

end


%% Plot PES

cmap = lines(num_pts);
hold on
for i = 1:4
    g{i} = plot(CC_dist,ENERGY_CC(i,:),'ko','markersize',5,'markerfacecolor',0.8*cmap(i,:),'color',0.8*cmap(i,:));
    plot(CC_dist,ENERGY_CC(i,:),'k-.','Linewidth',2,'color',0.8*cmap(i,:));

    %h{i} = plot(OH_dist,ENERGY_CRCC(i,:),'ko','markersize',5,'markerfacecolor',cmap(i,:),'color',cmap(i,:));
    %plot(OH_dist,ENERGY_CRCC(i,:),'k-','Linewidth',1.5,'color',cmap(i,:));
end
hold off
%ll1 = legend([h{1},h{2},h{3},h{4}],{'CR-CC(2,3)','S_1 CR-EOMCC(2,3)','S_2 CR-EOMCC(2,3)','S_3 CR-EOMCC(2,3)'}); set(ll1,'Location','NorthEast');
ll2 = legend([g{1},g{2},g{3},g{4}],{'CCSD','S_1 EOMCCSD','S_2 EOMCCSD','S_3 EOMCCSD'}); set(ll2,'Location','NorthEast');
xlabel('R_{OH} Distance / bohr')
ylabel('Energy / Eh')
set(gca,'FontSize',14,'Linewidth',2,'Box','off')
axis([-inf,inf,-76.6,-74])
grid on
