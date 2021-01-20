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
nact_h = 3; nact_p = 6; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
flag_act_scheme = 0;

sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                           nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-11;
ccopts.shift = 0;
flag_test = false;

% Full CCSDt(II): Ecorr = -0.220264751342 Eh / -0.616556843101 Eh
% CCSDt(II) in active scheme: Ecorr = -0.219121119129 Eh / -0.605291997421 Eh

%[cc_t,Ecorr] = uccsdt(sys,ccopts);
%[cc_t,Ecorr] = uccsd(sys,ccopts);
[cc_t,Ecorr] = uccsdt2(sys,ccopts,[],flag_test);
%[cc_t,Ecorr] = uccsdt2_full(sys,ccopts);

%% CC(t;3) correction

[HBar_t] = build_ucc_HBar( cc_t, sys, false);

lccopts.diis_size = 5;
lccopts.maxit = 500;
lccopts.tol = 1e-11;
lccopts.ground_shift = 0.0;

[cc_t,lccsd_resid] = luccsd(cc_t,HBar_t,sys,lccopts);

%
mm3_correction(cc_t,HBar_t,nact_h,nact_p,sys)


function [] = mm3_correction(cc_t,HBar_t,Nact_h,Nact_p,sys)

        %fprintf('\n')
        %fprintf('-------------------------------------------------------------------------------------\n')
        %fprintf('Performing correction for ground state -\n')
        
        Ecorr_ccsdt = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);
        Eccsdt = sys.Escf + Ecorr_ccsdt;
        %ed_result = -0.002414558899731;
                
        [deltaA,deltaB,deltaC,deltaD] = cct3_2(cc_t,HBar_t,Nact_h,Nact_p,sys);
        Ecorr_crcct3A = Ecorr_ccsdt + deltaA;
        Ecorr_crcct3B = Ecorr_ccsdt + deltaB;
        Ecorr_crcct3C = Ecorr_ccsdt + deltaC;
        Ecorr_crcct3D = Ecorr_ccsdt + deltaD;
        Ecrcct3A(1) = sys.Escf + Ecorr_crcct3A;
        Ecrcct3B(1) = sys.Escf + Ecorr_crcct3B;
        Ecrcct3C(1) = sys.Escf + Ecorr_crcct3C;
        Ecrcct3D(1) = sys.Escf + Ecorr_crcct3D;
        fprintf('\n')
        fprintf('CCSDt(II) = %4.12f Eh\n',Eccsdt)
        fprintf('CR-UCC(t,3)_A = %4.12f Eh     Ecorr_A = %4.12f Eh     Delta_A = %4.12f Eh\n',Ecrcct3A(1),Ecorr_crcct3A,deltaA)
        fprintf('CR-UCC(t,3)_B = %4.12f Eh     Ecorr_B = %4.12f Eh     Delta_B = %4.12f Eh\n',Ecrcct3B(1),Ecorr_crcct3B,deltaB)
        fprintf('CR-UCC(t,3)_C = %4.12f Eh     Ecorr_C = %4.12f Eh     Delta_C = %4.12f Eh\n',Ecrcct3C(1),Ecorr_crcct3C,deltaC)
        fprintf('CR-UCC(t,3)_D = %4.12f Eh     Ecorr_D = %4.12f Eh     Delta_D = %4.12f Eh\n',Ecrcct3D(1),Ecorr_crcct3D,deltaD)
        fprintf('-------------------------------------------------------------------------------------\n')

end