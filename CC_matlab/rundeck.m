clear all
clc
close all

addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests'));

%load h2o-631g-stretched
load chplus-1.0-olsen.mat

nfzc = 0; nfzv = 0; 
nact_h = 100; nact_p = 100; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
flag_act_scheme = 0;

sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                           nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);
                       
nroot = 10;     

ccopts.diis_size = 5;
ccopts.maxit = 200;
ccopts.tol = 1e-8;
ccopts.shift = 0;
                       
lccopts.maxit = 200;
lccopts.diis_size = 5;
lccopts.tol = 1e-7; 
lccopts.excited_shift = 0.0;
lccopts.ground_shift = 0.0;
lccopts.solver = 1;
lccopts.nroot = nroot;

eomopts.nroot = nroot;
eomopts.maxit = 200;
eomopts.tol = 1e-7;
eomopts.nvec_per_root = 1;
eomopts.max_nvec_per_root = 5;
eomopts.flag_verbose = 1;
eomopts.init_guess = 'cis';
eomopts.mult = 1;
eomopts.thresh_vec = 1e-3;
eomopts.solver = 2;
                      
run_crcc23(sys,nroot,ccopts,eomopts,lccopts)                     