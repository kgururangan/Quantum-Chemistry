clear all
clc
close all

format long

addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests'));

load h2o-631g-stretched
%load f2-2.0-pvdz.mat % NEEDS NFZC = 2

%% Set up problem

nfzc = 0; nfzv = 0; 
nact_h = 2; nact_p = 4; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
flag_act_scheme = 0;

sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                           nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);
                       
PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;


                       
%% Full CCSDT

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-10;
ccopts.shift = 0;

[cc_t,Ecorr_uccsdt] = uccsdt(sys,ccopts);

%% t2B - projection 1 (PPHH) [DONE]
clc

out_proj = 'ABIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnIf),t3b(AfBmnJ)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfJ),t3b(AfBInm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnJf),t3c(AfBInm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnIf),t3c(AfBmnJ)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(Anef),t3b(efBInJ)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(Anef),t3c(efBInJ)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nBfe),t3b(AfeInJ)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(Bnef),t3c(AfeInJ)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(AeBImJ)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(AeBImJ)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_1] = update_t2b_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_PPHH = %4.30f\n',get_error(X2B(PA,PB,HA,HB),X2B_1(PA,PB,HA,HB)))

%% t2B - projection 2 (PpHH) [DONE]
clc

out_proj = 'AbIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnIf),t3b(AfbmnJ)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfJ),t3b(AfbInm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnJf),t3c(AfbInm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnIf),t3c(AfbmnJ)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(Anef),t3b(efbInJ)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(Anef),t3c(efbInJ)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nbfe),t3b(AfeInJ)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(bnef),t3c(AfeInJ)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(AebImJ)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(AebImJ)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_2] = update_t2b_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_PpHH = %4.15f\n',get_error(X2B(PA,pB,HA,HB),X2B_2(PA,pB,HA,HB)))

%% t2B - projection 3 (pPHH) [DONE]
clc

out_proj = 'aBIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnIf),t3b(afBmnJ)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfJ),t3b(afBInm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnJf),t3c(afBInm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnIf),t3c(afBmnJ)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(anef),t3b(efBInJ)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(anef),t3c(efBInJ)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nBfe),t3b(afeInJ)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(Bnef),t3c(afeInJ)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(aeBImJ)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(aeBImJ)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_3] = update_t2b_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_pPHH = %4.15f\n',get_error(X2B(pA,PB,HA,HB),X2B_3(pA,PB,HA,HB)))

%% t2B - projection 4 (PPhH) [DONE]
clc

out_proj = 'ABiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnif),t3b(AfBmnJ)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfJ),t3b(AfBinm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnJf),t3c(AfBinm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnif),t3c(AfBmnJ)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(Anef),t3b(efBinJ)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(Anef),t3c(efBinJ)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nBfe),t3b(AfeinJ)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(Bnef),t3c(AfeinJ)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(AeBimJ)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(AeBimJ)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_4] = update_t2b_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_PPhH = %4.15f\n',get_error(X2B(PA,PB,hA,HB),X2B_4(PA,PB,hA,HB)))

%% t2B - projection 5 (PPHh) [DONE]
clc

out_proj = 'ABIj';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnIf),t3b(AfBmnj)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfj),t3b(AfBInm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnjf),t3c(AfBInm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnIf),t3c(AfBmnj)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(Anef),t3b(efBInj)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(Anef),t3c(efBInj)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nBfe),t3b(AfeInj)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(Bnef),t3c(AfeInj)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(AeBImj)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(AeBImj)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_5] = update_t2b_proj5(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_PPHh = %4.15f\n',get_error(X2B(PA,PB,HA,hB),X2B_5(PA,PB,HA,hB)))

%% t2B - projection 6 (PphH) [DONE]
clc

out_proj = 'AbiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnif),t3b(AfbmnJ)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfJ),t3b(Afbinm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnJf),t3c(Afbinm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnif),t3c(AfbmnJ)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(Anef),t3b(efbinJ)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(Anef),t3c(efbinJ)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nbfe),t3b(AfeinJ)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(bnef),t3c(AfeinJ)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(AebimJ)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(AebimJ)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_6] = update_t2b_proj6(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_PphH = %4.15f\n',get_error(X2B(PA,pB,hA,HB),X2B_6(PA,pB,hA,HB)))

%% t2B - projection 7 (pPhH) [DONE]
clc

out_proj = 'aBiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnif),t3b(afBmnJ)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfJ),t3b(afBinm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnJf),t3c(afBinm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnif),t3c(afBmnJ)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(anef),t3b(efBinJ)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(anef),t3c(efBinJ)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nBfe),t3b(afeinJ)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(Bnef),t3c(afeinJ)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(aeBimJ)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(aeBimJ)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_7] = update_t2b_proj7(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_pPhH = %4.15f\n',get_error(X2B(pA,PB,hA,HB),X2B_7(pA,PB,hA,HB)))

%% t2B - projection 8 (PpHh) [DONE]
clc

out_proj = 'AbIj';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnIf),t3b(Afbmnj)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfj),t3b(AfbInm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnjf),t3c(AfbInm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnIf),t3c(Afbmnj)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(Anef),t3b(efbInj)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(Anef),t3c(efbInj)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nbfe),t3b(AfeInj)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(bnef),t3c(AfeInj)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(AebImj)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(AebImj)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_8] = update_t2b_proj8(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_PpHh = %4.15f\n',get_error(X2B(PA,pB,HA,hB),X2B_8(PA,pB,HA,hB)))

%% t2B - projection 9 (pPHh) [DONE]
clc

out_proj = 'aBIj';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnIf),t3b(afBmnj)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfj),t3b(afBInm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnjf),t3c(afBInm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnIf),t3c(afBmnj)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(anef),t3b(efBInj)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(anef),t3c(efBInj)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nBfe),t3b(afeInj)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(Bnef),t3c(afeInj)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(aeBImj)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(aeBImj)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_9] = update_t2b_proj9(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_pPHh = %4.15f\n',get_error(X2B(pA,PB,HA,hB),X2B_9(pA,PB,HA,hB)))

%% t2B - projection 10 (PPhh) [DONE]
clc

out_proj = 'ABij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnif),t3b(AfBmnj)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfj),t3b(AfBinm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnjf),t3c(AfBinm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnif),t3c(AfBmnj)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(Anef),t3b(efBinj)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(Anef),t3c(efBinj)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nBfe),t3b(Afeinj)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(Bnef),t3c(Afeinj)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(AeBimj)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(AeBimj)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_10] = update_t2b_proj10(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_PPhh = %4.15f\n',get_error(X2B(PA,PB,hA,hB),X2B_10(PA,PB,hA,hB)))


%% t2B - projection 11 (ppHH) [DONE]
clc

out_proj = 'abIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnIf),t3b(afbmnJ)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfJ),t3b(afbInm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnJf),t3c(afbInm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnIf),t3c(afbmnJ)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(anef),t3b(efbInJ)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(anef),t3c(efbInJ)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nbfe),t3b(afeInJ)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(bnef),t3c(afeInJ)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(aebImJ)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(aebImJ)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_11] = update_t2b_proj11(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_ppHH = %4.15f\n',get_error(X2B(pA,pB,HA,HB),X2B_11(pA,pB,HA,HB)))

%% t2B - projection 12 (Pphh) [DONE]
clc

out_proj = 'Abij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnif),t3b(Afbmnj)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfj),t3b(Afbinm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnjf),t3c(Afbinm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnif),t3c(Afbmnj)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(Anef),t3b(efbinj)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(Anef),t3c(efbinj)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nbfe),t3b(Afeinj)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(bnef),t3c(Afeinj)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(Aebimj)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(Aebimj)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_12] = update_t2b_proj12(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_Pphh = %4.15f\n',get_error(X2B(PA,pB,hA,hB),X2B_12(PA,pB,hA,hB)))

%% t2B - projection 13 (pPhh) [DONE]
clc

out_proj = 'aBij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnif),t3b(afBmnj)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfj),t3b(afBinm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnjf),t3c(afBinm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnif),t3c(afBmnj)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(anef),t3b(efBinj)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(anef),t3c(efBinj)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nBfe),t3b(afeinj)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(Bnef),t3c(afeinj)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(aeBimj)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(aeBimj)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_13] = update_t2b_proj13(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_pPhh = %4.15f\n',get_error(X2B(pA,PB,hA,hB),X2B_13(pA,PB,hA,hB)))


%% t2B - projection 14 (pphH) [DONE]
clc

out_proj = 'abiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnif),t3b(afbmnJ)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfJ),t3b(afbinm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnJf),t3c(afbinm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnif),t3c(afbmnJ)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(anef),t3b(efbinJ)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(anef),t3c(efbinJ)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nbfe),t3b(afeinJ)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(bnef),t3c(afeinJ)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(aebimJ)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(aebimJ)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_14] = update_t2b_proj14(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_pphH = %4.15f\n',get_error(X2B(pA,pB,hA,HB),X2B_14(pA,pB,hA,HB)))

%% t2B - projection 15 (ppHh) [DONE]
clc

out_proj = 'abIj';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnIf),t3b(afbmnj)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfj),t3b(afbInm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnjf),t3c(afbInm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnIf),t3c(afbmnj)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(anef),t3b(efbInj)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(anef),t3c(efbInj)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nbfe),t3b(afeInj)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(bnef),t3c(afeInj)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(aebImj)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(aebImj)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_15] = update_t2b_proj15(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_ppHh = %4.15f\n',get_error(X2B(pA,pB,HA,hB),X2B_15(pA,pB,HA,hB)))

%% t2B - projection 16 (pphh) [DONE]
clc

out_proj = 'abij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = -0.5*einsum_kg(H2A.ooov,t3b,'mnif,afbmnj->abij');
d1 = {'-h2A(mnif),t3b(afbmnj)'};
fun(d1,'D1')

%TR_D2 = -einsum_kg(H2B.oovo,t3b,'nmfj,afbinm->abij');
d2 = {'-h2B(nmfj),t3b(afbinm)'};
fun(d2,'D2')

%TR_D3 = -0.5*einsum_kg(H2C.ooov,t3c,'mnjf,afbinm->abij');
d3 = {'-h2C(mnjf),t3c(afbinm)'};
fun(d3,'D3')

%TR_D4 = -einsum_kg(H2B.ooov,t3c,'mnif,afbmnj->abij');
d4 = {'-h2B(mnif),t3c(afbmnj)'};
fun(d4,'D4')

%TR_D5 = +0.5*einsum_kg(H2A.vovv,t3b,'anef,efbinj->abij');
d5 = {'h2A(anef),t3b(efbinj)'};
fun(d5,'D5')

%TR_D6 = +einsum_kg(H2B.vovv,t3c,'anef,efbinj->abij');
d6 = {'h2B(anef),t3c(efbinj)'};
fun(d6,'D6')

%TR_D7 = +einsum_kg(H2B.ovvv,t3b,'nbfe,afeinj->abij');
d7 = {'h2B(nbfe),t3b(afeinj)'};
fun(d7,'D7')

%TR_D8 = +0.5*einsum_kg(H2C.vovv,t3c,'bnef,afeinj->abij');
d8 = {'h2C(bnef),t3c(afeinj)'};
fun(d8,'D8')

%TR_D9 = +einsum_kg(H1A.ov,t3b,'me,aebimj->abij');
d9 = {'h1A(me),t3b(aebimj)'};
fun(d9,'D9')

%TR_D10 = +einsum_kg(H1B.ov,t3c,'me,aebimj->abij');
d10 = {'h1B(me),t3c(aebimj)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2B] = build_t2b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2B_16] = update_t2b_proj16(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2B_pphh = %4.15f\n',get_error(X2B(pA,pB,hA,hB),X2B_16(pA,pB,hA,hB)))

%%

function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end

