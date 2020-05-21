clear all
clc
close all

format long

addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/Davidson_matlab_test'));
addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/GDav'));
addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/Bchdav'));

addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2o-pvdz/'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2o-631g/'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2-pvdz/'));


%load h2o-pvdz-integrals.mat
%load h2-pvdz-integrals.mat
%load h2o-631g-stretched-integrals
load h2o-631g-eom-ccpy

sys_ucc = build_system_ucc(e1int,e2int,Nocc_a,Nocc_b);
sys_cc = build_system_cc(e1int,e2int,Nocc_a+Nocc_b);


% load h2o-pvdz-ccpy.mat
%sys_cc = build_system(VM,FM,occ,unocc,2);

%% UCCSD

% update: 4/6/20 - 
% there is still something wrong with UCC - errors on order of 10^-6... but accurate for the most part...
% problem seems to lie primarily in T2B with some errors in T2A and T2C

ccopts.diis_size = 3;
ccopts.maxit = 100;
ccopts.tol = 1e-8;
[t1a,t1b,t2a,t2b,t2c,Ecorr_ucc] = uccsd(sys_ucc,ccopts);
[T1] = convert_spinint_to_spinorb({t1a,t1b},sys_ucc);
[T2] = convert_spinint_to_spinorb({t2a,t2b,t2c},sys_ucc);

%% Checking UCCSD results with Jun's

E1A_exact = 9.944089519854425E-011;
E1B_exact = 9.944091634197451E-011;
E2A_exact =    -1.118149855483151E-002;
E2C_exact =    -1.118149857437544E-002;
E2B_exact =    -0.193074895116479;
E1A1A_exact =   -7.721105018602925E-005;
E1A1B_exact =   5.956425520525440E-004;
E1B1B_exact =   -7.721105255510347E-005;
Ecorr_exact = E1A_exact + E1B_exact + E2A_exact + E2B_exact + E2C_exact + E1A1A_exact + E1B1B_exact + E1A1B_exact;

E1A = einsum(sys_ucc.fa_ov,t1a,'me,em->');
E1B = einsum(sys_ucc.fb_ov,t1b,'me,em->');
E2A = 0.25*einsum(sys_ucc.vA_oovv,t2a,'mnef,efmn->');
E2C = 0.25*einsum(sys_ucc.vC_oovv,t2c,'mnef,efmn->');
E2B = einsum(sys_ucc.vB_oovv,t2b,'mnef,efmn->');
E1A1A = 0.5*einsum(einsum(sys_ucc.vA_oovv,t1a,'mnef,fn->me'),t1a,'me,em->');
E1B1B = 0.5*einsum(einsum(sys_ucc.vC_oovv,t1b,'mnef,fn->me'),t1b,'me,em->');
E1A1B = einsum(einsum(sys_ucc.vB_oovv,t1a,'mnef,em->nf'),t1b,'nf,fn->');

Ecorr = E1A + E1B + E2A + E2B + E2C + E1A1A + E1B1B + E1A1B;

fprintf('E1A error = %4.12f\n',E1A_exact - E1A)
fprintf('E1B error = %4.12f\n',E1B_exact - E1B)
fprintf('E2A error = %4.12f\n',E2A_exact - E2A)
fprintf('E2C error = %4.12f\n',E2C_exact - E2C)
fprintf('E2B error = %4.12f\n',E2B_exact - E2B)
fprintf('E1A1A error = %4.12f\n',E1A1A_exact - E1A1A)
fprintf('E1B1B error = %4.12f\n',E1B1B_exact - E1B1B)
fprintf('E1A1B error = %4.12f\n',E1A1B_exact - E1A1B)
fprintf('Ecorr error = %4.12f\n',Ecorr_exact - Ecorr)

        
%% CCSD

ccopts.diis_size = 3;
ccopts.maxit = 100;
ccopts.tol = 1e-9;

[t1,t2,Ecorr_cc] = ccsd(sys_cc,ccopts);


%% Build HBar

[HBar] = build_HBar(t1,t2,sys_cc); 
[D, Dia, Dijab] = HBar_CCSD_diagonal(HBar, t1, t2, sys_cc);
[Nocc,Nunocc] = size(HBar{1}{1,2}); Nov = Nocc*Nunocc; Hbar_dim = Nov + Nov^2;

%% Left-CCSD ground state

ccopts.diis_size = 5;
ccopts.maxit = 200;
ccopts.tol = 1e-9;

[lambda1,lambda2,lcc_resid] = lccsd(0.0,[],t1,t2,HBar,sys_cc,ccopts);


%% EOM-CCSD

eomopts.nroot = 15;
eomopts.maxit = 1000;
eomopts.tol = 1e-7;
eomopts.nvec_per_root = 1;
eomopts.max_nvec_per_root = 10;
eomopts.flag_verbose = 1;
eomopts.init_guess = 'cis';
eomopts.thresh_vec = 1e-7; % normally set to the convergence threshold of eom

% it's possible there's an error in EOMCCSD sigma function. R vector for
% first singlet state of H2 are close to Jun's but not exactly identical...
[Rvec, omega, r0, eom_residual] = eomccsd(HBar,t1,t2,sys_cc,eomopts);

% Using iterative diagonalization solvers of Yunkai Zhou
% these also encounter subspace collapse and 0 eigenvalue problems!

% eomopts.tol = 1e-4;
% eomopts.Matsymm = 1;
% eomopts.vmax = max(2*nroot, 20);
% eomopts.itmax = 100;
% eomopts.nkeep = nroot;
% eomopts.v0 = rand(Nocc*Nunocc+Nocc^2*Nunocc^2,1);
% eomopts.displ = 2;
% eomopts.les_solver = 31; % 71 jacobi-davidson minres does not work!

% HRvec = @(x) HR_matvec(x,HBar);
% [eval, V, nconv, history] = gdav(HRvec, Nocc*Nunocc+Nocc^2*Nunocc^2, eomopts.nroot, 'SM',eomopts);
% [eval, V, nconv, history] = bchdav('HR_matvec', Nocc*Nunocc+Nocc^2*Nunocc^2, eomopts.nroot);

%% Left-EOMCCSD 

ccopts.maxit = 100;
ccopts.diis_size = 3;
ccopts.tol = 1e-6; % can only converge Left-EOM to the same tolerance as we converged right EOM

[Lvec,eomlcc_resid] = lefteomccsd(omega,Rvec,HBar,t1,t2,sys_cc,ccopts);



%%







