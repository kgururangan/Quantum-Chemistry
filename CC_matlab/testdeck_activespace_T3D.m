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

%% t3D - projection 1 (PPPHHH) [DONE]
clc

out_proj = 'ABCIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1B.oo,t3d,'mk,abcijm->abcijk');
d1 = {'-h1B(mK),t3d(ABCIJm)'};
fun(d1,'D1')

%D2 = einsum_kg(H1B.vv,t3d,'ce,abeijk->abcijk');
d2 = {'h1B(Ce),t3d(ABeIJK)'};
fun(d2,'D2')

%D3 = 0.5*einsum_kg(H2C.oooo,t3d,'mnij,abcmnk->abcijk');
d3 = {'h2C(mnIJ),t3d(ABCmnK)'};
fun(d3,'D3')

%D4 = 0.5*einsum_kg(H2C.vvvv,t3d,'abef,efcijk->abcijk');
d4 = {'h2C(ABef),t3d(efCIJK)'};
fun(d4,'D4')

%D5 = einsum_kg(H2B.ovvo,t3c,'maei,ebcmjk->abcijk');
d5 = {'h2B(mAeI),t3c(eBCmJK)'};
fun(d5,'D5')

%D6 = einsum_kg(H2C.voov,t3d,'amie,ebcmjk->abcijk');
d6 = {'h2C(AmIe),t3d(eBCmJK)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3D,VT3D] = build_t3d(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3D_1] = update_t3d_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3D_PPPHHH = %4.15f\n',get_error(X3D(PB,PB,PB,HB,HB,HB),X3D_1))

%% t3D - projection 2 (PPpHHH) [DONE]
clc

out_proj = 'ABcIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1B.oo,t3d,'mk,abcijm->abcijk');
d1 = {'-h1B(mK),t3d(ABcIJm)'};
fun(d1,'D1')

%D2 = einsum_kg(H1B.vv,t3d,'ce,abeijk->abcijk');
d2 = {'h1B(ce),t3d(ABeIJK)',};
fun(d2,'D2')
d3 = {'-h1B(Be),t3d(AceIJK)'};
fun(d3,'D3')

%D3 = 0.5*einsum_kg(H2C.oooo,t3d,'mnij,abcmnk->abcijk');
d4 = {'h2C(mnIJ),t3d(ABcmnK)'};
fun(d4,'D4')

%D4 = 0.5*einsum_kg(H2C.vvvv,t3d,'abef,efcijk->abcijk');
d5 = {'h2C(ABef),t3d(efcIJK)'};
fun(d5,'D5')
d6 = {'-h2C(Acef),t3d(efBIJK)'};
fun(d6,'D6')

%D5 = einsum_kg(H2B.ovvo,t3c,'maei,ebcmjk->abcijk');
d7 = {'h2B(mAeI),t3c(eBcmJK)'};
fun(d7,'D7')
d8 = {'-h2B(mceI),t3c(eBAmJK)'};
fun(d8,'D8')

%D6 = einsum_kg(H2C.voov,t3d,'amie,ebcmjk->abcijk');
d9 = {'h2C(AmIe),t3d(eBcmJK)'};
fun(d9,'D9')
d10 = {'-h2C(cmIe),t3d(eBAmJK)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3D,VT3D] = build_t3d(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3D_2] = update_t3d_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3D_PPpHHH = %4.15f\n',get_error(X3D(PB,PB,pB,HB,HB,HB),X3D_2))

%% t3D - projection 3 (PPPhHH) [DONE]
clc

out_proj = 'ABCiJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1B.oo,t3d,'mk,abcijm->abcijk');
d1 = {'-h1B(mK),t3d(ABCiJm)'};
fun(d1,'D1')
d2 = {'h1B(mi),t3d(ABCKJm)'};
fun(d2,'D2')

%D2 = einsum_kg(H1B.vv,t3d,'ce,abeijk->abcijk');
d3 = {'h1B(Ce),t3d(ABeiJK)'};
fun(d3,'D3')

%D3 = 0.5*einsum_kg(H2C.oooo,t3d,'mnij,abcmnk->abcijk');
d4 = {'h2C(mniJ),t3d(ABCmnK)'};
fun(d4,'D4')
d5 = {'-h2C(mnKJ),t3d(ABCmni)'};
fun(d5,'D5')

%D4 = 0.5*einsum_kg(H2C.vvvv,t3d,'abef,efcijk->abcijk');
d6 = {'h2C(ABef),t3d(efCiJK)'};
fun(d6,'D6')

%D5 = einsum_kg(H2B.ovvo,t3c,'maei,ebcmjk->abcijk');
d7 = {'h2B(mAei),t3c(eBCmJK)'};
fun(d7,'D7')
d8 = {'-h2B(mAeJ),t3c(eBCmiK)'};
fun(d8,'D8')

%D6 = einsum_kg(H2C.voov,t3d,'amie,ebcmjk->abcijk');
d9 = {'h2C(Amie),t3d(eBCmJK)'};
fun(d9,'D9')
d10 = {'-h2C(AmJe),t3d(eBCmiK)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3D,VT3D] = build_t3d(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3D_3] = update_t3d_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3D_PPPhHH = %4.15f\n',get_error(X3D(PB,PB,PB,hB,HB,HB),X3D_3))

%% t3D - projection 4 (PPphHH) [DONE]
clc

out_proj = 'ABciJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1B.oo,t3d,'mk,abcijm->abcijk');
d1 = {'-h1B(mK),t3d(ABciJm)'};
fun(d1,'D1')
d2 = {'h1B(mi),t3d(ABcKJm)'};
fun(d2,'D2')

%D2 = einsum_kg(H1B.vv,t3d,'ce,abeijk->abcijk');
d3 = {'h1B(ce),t3d(ABeiJK)'};
fun(d3,'D3')
d4 = {'-h1B(Be),t3d(AceiJK)'};
fun(d4,'D4')

%D3 = 0.5*einsum_kg(H2C.oooo,t3d,'mnij,abcmnk->abcijk');
d5 = {'h2C(mniJ),t3d(ABcmnK)'};
fun(d5,'D5')
d6 = {'-h2C(mnKJ),t3d(ABcmni)'};
fun(d6,'D6')

%D4 = 0.5*einsum_kg(H2C.vvvv,t3d,'abef,efcijk->abcijk');
d7 = {'h2C(ABef),t3d(efciJK)'};
fun(d7,'D7')
d8 = {'-h2C(Acef),t3d(efBiJK)'};
fun(d8,'D8')

%D5 = einsum_kg(H2B.ovvo,t3c,'maei,ebcmjk->abcijk');
d9 = {'h2B(mAei),t3c(eBcmJK)'};
fun(d9,'D9')
d10 = {'-h2B(mcei),t3c(eBAmJK)'};
fun(d10,'D10')
d11 = {'-h2B(mAeJ),t3c(eBcmiK)'};
fun(d11,'D11')
d12 = {'h2B(mceJ),t3c(eBAmiK)'};
fun(d12,'D12')

%D6 = einsum_kg(H2C.voov,t3d,'amie,ebcmjk->abcijk');
d13 = {'h2C(Amie),t3d(eBcmJK)'};
fun(d13,'D13')
d14 = {'-h2C(cmie),t3d(eBAmJK)'};
fun(d14,'D14')
d15 = {'-h2C(AmJe),t3d(eBcmiK)'};
fun(d15,'D15')
d16 = {'h2C(cmJe),t3d(eBAmiK)'};
fun(d16,'D16')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3D,VT3D] = build_t3d(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3D_4] = update_t3d_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3D_PPphHH = %4.15f\n',get_error(X3D(PB,PB,pB,hB,HB,HB),X3D_4))


%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end