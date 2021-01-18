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

%% t1A - projection 1 (PH)
clc

out_proj = 'AI';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = +0.25*einsum_kg(sys.vA_oovv,t3a,'mnef,aefimn->ai');
d1 = {'v2A(mnef),t3a(AefImn)'};
fun(d1,'D1')

% TR_D2 = +einsum_kg(sys.vB_oovv,t3b,'mnef,aefimn->ai'); 
d2 = {'v2B(mnef),t3b(AefImn)'};
fun(d2,'D2')
 
% TR_D3 = +0.25*einsum_kg(sys.vC_oovv,t3c,'mnef,aefimn->ai');
d3 = {'v2C(mnef),t3c(AefImn)'};
fun(d3,'D3')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X1A] = build_t1a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X1A_1] = update_t1a_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X1A_PH = %4.15f\n',get_error(X1A(PA,HA),X1A_1(PA,HA)))

%% t1A - projection 2 (pH)
clc

out_proj = 'aI';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = +0.25*einsum_kg(sys.vA_oovv,t3a,'mnef,aefimn->ai');
d1 = {'v2A(mnef),t3a(aefImn)'};
fun(d1,'D1')

% TR_D2 = +einsum_kg(sys.vB_oovv,t3b,'mnef,aefimn->ai'); 
d2 = {'v2B(mnef),t3b(aefImn)'};
fun(d2,'D2')
 
% TR_D3 = +0.25*einsum_kg(sys.vC_oovv,t3c,'mnef,aefimn->ai');
d3 = {'v2C(mnef),t3c(aefImn)'};
fun(d3,'D3')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X1A] = build_t1a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X1A_2] = update_t1a_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X1A_pH = %4.15f\n',get_error(X1A(pA,HA),X1A_2(pA,HA)))

%% t1A - projection 3 (Ph)
clc

out_proj = 'Ai';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = +0.25*einsum_kg(sys.vA_oovv,t3a,'mnef,aefimn->ai');
d1 = {'v2A(mnef),t3a(Aefimn)'};
fun(d1,'D1')

% TR_D2 = +einsum_kg(sys.vB_oovv,t3b,'mnef,aefimn->ai'); 
d2 = {'v2B(mnef),t3b(Aefimn)'};
fun(d2,'D2')
 
% TR_D3 = +0.25*einsum_kg(sys.vC_oovv,t3c,'mnef,aefimn->ai');
d3 = {'v2C(mnef),t3c(Aefimn)'};
fun(d3,'D3')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X1A] = build_t1a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X1A_3] = update_t1a_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X1A_Ph = %4.15f\n',get_error(X1A(PA,hA),X1A_3(PA,hA)))

%% t1A - projection 4 (ph)
clc

out_proj = 'ai';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = +0.25*einsum_kg(sys.vA_oovv,t3a,'mnef,aefimn->ai');
d1 = {'v2A(mnef),t3a(aefimn)'};
fun(d1,'D1')

% TR_D2 = +einsum_kg(sys.vB_oovv,t3b,'mnef,aefimn->ai'); 
d2 = {'v2B(mnef),t3b(aefimn)'};
fun(d2,'D2')
 
% TR_D3 = +0.25*einsum_kg(sys.vC_oovv,t3c,'mnef,aefimn->ai');
d3 = {'v2C(mnef),t3c(aefimn)'};
fun(d3,'D3')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X1A] = build_t1a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X1A_4] = update_t1a_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X1A_ph = %4.15f\n',get_error(X1A(pA,hA),X1A_4(pA,hA)))

%% t1B - projection 1 (PH)
clc

out_proj = 'AI';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = +0.25*einsum_kg(sys.vC_oovv,t3d,'mnef,aefimn->ai');
d1 = {'v2C(mnef),t3d(AefImn)'};
fun(d1,'D1')

% TR_D2 = +einsum_kg(sys.vB_oovv,t3c,'mnef,efamni->ai'); 
d2 = {'v2B(mnef),t3c(efAmnI)'};
fun(d2,'D2')
 
% TR_D3 = +0.25*einsum_kg(sys.vA_oovv,t3b,'mnef,efamni->ai');
d3 = {'v2A(mnef),t3b(efAmnI)'};
fun(d3,'D3')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X1B] = build_t1b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X1B_1] = update_t1b_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X1B_PH = %4.15f\n',get_error(X1B(PB,HB),X1B_1(PB,HB)))

%% t1B - projection 2 (pH)
clc

out_proj = 'aI';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = +0.25*einsum_kg(sys.vC_oovv,t3d,'mnef,aefimn->ai');
d1 = {'v2C(mnef),t3d(aefImn)'};
fun(d1,'D1')

% TR_D2 = +einsum_kg(sys.vB_oovv,t3c,'mnef,efamni->ai'); 
d2 = {'v2B(mnef),t3c(efamnI)'};
fun(d2,'D2')
 
% TR_D3 = +0.25*einsum_kg(sys.vA_oovv,t3b,'mnef,efamni->ai');
d3 = {'v2A(mnef),t3b(efamnI)'};
fun(d3,'D3')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X1B] = build_t1b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X1B_2] = update_t1b_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X1B_pH = %4.15f\n',get_error(X1B(pB,HB),X1B_2(pB,HB)))

%% t1B - projection 3 (Ph)
clc

out_proj = 'Ai';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = +0.25*einsum_kg(sys.vC_oovv,t3d,'mnef,aefimn->ai');
d1 = {'v2C(mnef),t3d(Aefimn)'};
fun(d1,'D1')

% TR_D2 = +einsum_kg(sys.vB_oovv,t3c,'mnef,efamni->ai'); 
d2 = {'v2B(mnef),t3c(efAmni)'};
fun(d2,'D2')
 
% TR_D3 = +0.25*einsum_kg(sys.vA_oovv,t3b,'mnef,efamni->ai');
d3 = {'v2A(mnef),t3b(efAmni)'};
fun(d3,'D3')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X1B] = build_t1b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X1B_3] = update_t1b_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X1B_Ph = %4.15f\n',get_error(X1B(PB,hB),X1B_3(PB,hB)))

%% t1B - projection 1 (PH)
clc

out_proj = 'ai';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = +0.25*einsum_kg(sys.vC_oovv,t3d,'mnef,aefimn->ai');
d1 = {'v2C(mnef),t3d(aefimn)'};
fun(d1,'D1')

% TR_D2 = +einsum_kg(sys.vB_oovv,t3c,'mnef,efamni->ai'); 
d2 = {'v2B(mnef),t3c(efamni)'};
fun(d2,'D2')
 
% TR_D3 = +0.25*einsum_kg(sys.vA_oovv,t3b,'mnef,efamni->ai');
d3 = {'v2A(mnef),t3b(efamni)'};
fun(d3,'D3')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X1B] = build_t1b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X1B_4] = update_t1b_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X1B_ph = %4.15f\n',get_error(X1B(pB,hB),X1B_4(pB,hB)))

%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end
