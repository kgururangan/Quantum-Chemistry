clear all
clc
close all

format long

addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/CAD_triples_correction/tests/cct3_tests/act-2-2'))

%load emiliano-h2o-1.0-test-ints.mat

load h2o-631g-stretched
%load f2-2.0-pvdz.mat % NEEDS NFZC = 2

nfzc = 0; nfzv = 0; 
nact_h = 2; nact_p = 6; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
flag_act_scheme = 0;

sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                           nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);

ccopts.diis_size = 5;
ccopts.maxit = 50;
ccopts.tol = 1e-11;
ccopts.shift = 0;
flag_test = false;

%[cc_t,Ecorr] = uccsdt(sys,ccopts);
%[cc_t,Ecorr] = uccsd(sys,ccopts);
[cc_t,Ecorr] = uccsdt2(sys,ccopts,[],flag_test);
%[cc_t,Ecorr] = uccsdt2_full(sys,ccopts);

