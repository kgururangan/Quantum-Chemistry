clear all
clc
close all

format long

addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
%addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/Davidson_matlab_test'));
%addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/GDav'));
%addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/Bchdav'));

addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2o-pvdz/'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2o-631g/'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2-pvdz/'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/'));

%%
workpath = '/Users/karthik/Desktop/CC_matlab_tests/h2o-631g/stretched/';
[e1int, e2int, Vnuc, Norb] = load_integrals(workpath);
Nocc_a = 5;
Nocc_b = 5;
Nelec = 10;

%%
%load h2o-pvdz.mat
%load h2-pvdz-integrals.mat
load h2o-631g-stretched
%load h2o-631g
%load rectangle-pvdz-d2h

nfzc = 0; nfzv = 0; nact_h = 200; nact_p = 200;

sys_ucc = build_system_ucc(e1int,e2int,Nocc_a,Nocc_b);
sys_cc = build_system_cc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,nfzv,nact_h,nact_p);

%% UCCSD

ccopts.diis_size = 3;
ccopts.maxit = 100;
ccopts.tol = 1e-10;
[t1a,t1b,t2a,t2b,t2c,Ecorr_ucc] = uccsd(sys_ucc,ccopts);
[T1] = convert_spinint_to_spinorb({t1a,t1b},sys_ucc);
[T2] = convert_spinint_to_spinorb({t2a,t2b,t2c},sys_ucc);

% % Checking UCCSD results with Jun's

E1A_exact = 9.944089519854425E-011;
E1B_exact = 9.944091634197451E-011;
E2A_exact =    -1.118149855483151E-002;
E2C_exact =    -1.118149857437544E-002;
E2B_exact =    -0.193074895116479;
E1A1A_exact =   -7.721105018602925E-005;
E1A1B_exact =   5.956425520525440E-004;
E1B1B_exact =   -7.721105255510347E-005;
Ecorr_exact = E1A_exact + E1B_exact + E2A_exact + E2B_exact + E2C_exact + E1A1A_exact + E1B1B_exact + E1A1B_exact;

E1A = einsum_kg(sys_ucc.fa_ov,t1a,'me,em->');
E1B = einsum_kg(sys_ucc.fb_ov,t1b,'me,em->');
E2A = 0.25*einsum_kg(sys_ucc.vA_oovv,t2a,'mnef,efmn->');
E2C = 0.25*einsum_kg(sys_ucc.vC_oovv,t2c,'mnef,efmn->');
E2B = einsum_kg(sys_ucc.vB_oovv,t2b,'mnef,efmn->');
E1A1A = 0.5*einsum_kg(einsum_kg(sys_ucc.vA_oovv,t1a,'mnef,fn->me'),t1a,'me,em->');
E1B1B = 0.5*einsum_kg(einsum_kg(sys_ucc.vC_oovv,t1b,'mnef,fn->me'),t1b,'me,em->');
E1A1B = einsum_kg(einsum_kg(sys_ucc.vB_oovv,t1a,'mnef,em->nf'),t1b,'nf,fn->');

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
ccopts.tol = 1e-8;

[t1,t2,Ecorr_ccsd] = ccsd(sys_cc,ccopts);

%% CCSDT

% CCSDT errors on order of 10^-5... be VERY careful checking all antisymmetrizers with einsum!
[t1,t2,t3,Ecorr_ccsdt] = ccsdt(sys_cc,ccopts);
%[t1,t2,t3,Ecorr_ccsdt3] = ccsdt3(sys_cc,ccopts);


%% Build HBar

%[HBar] = build_HBar(t1,t2,sys_cc); 
[HBar] = build_HBar_debug(t1,t2,sys_cc);

%% Left-CCSD 

ccopts.diis_size = 5;
ccopts.maxit = 200;
ccopts.tol = 1e-9;

[lambda1,lambda2,lcc_resid] = lccsd(t1,t2,HBar,sys_cc,ccopts);

%% CR-CC(2,3)

%[crcc23A, crcc23B, crcc23C, crcc23D] = crcc23(t1,t2,lambda1,lambda2,HBar,sys_cc);
[Ecorr_crcc23A,Ecorr_crcc23B,Ecorr_crcc23C,Ecorr_crcc23D] = crcc23_opt(t1,t2,lambda1,lambda2,HBar,sys_cc);

%% EOM-CCSD

eomopts.nroot = 3;
eomopts.maxit = 100;
eomopts.tol = 1e-6;
eomopts.nvec_per_root = 1;
eomopts.max_nvec_per_root = 5;
eomopts.flag_verbose = 1;
eomopts.init_guess = 'cis';
eomopts.thresh_vec = 10*eomopts.tol;

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

lccopts.maxit = 100;
lccopts.diis_size = 3;
lccopts.tol = eomopts.tol; % can only converge Left-EOM to the same tolerance as we converged right EOM

[Lvec,omega_lcc,eom_lcc_resid] = lefteomccsd(omega,Rvec,HBar,t1,t2,sys_cc,lccopts);

