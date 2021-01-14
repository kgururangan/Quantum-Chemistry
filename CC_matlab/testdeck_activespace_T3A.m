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

%% t3A - projection 1 (PPPHHH) [DONE]
clc

out_proj = 'ABCIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
d1 = {'-h1A(mK),t3a(ABCIJm)'};
fun(d1,'D1')

% D2 = einsum_kg(H1A.vv,t3a,'ce,abeijk->abcijk');
d2 = {'h1A(Ce),t3a(ABeIJK)'};
fun(d2,'D2')

% D3 = 0.5*einsum_kg(H2A.oooo,t3a,'mnij,abcmnk->abcijk');
d3 = {'h2A(mnIJ),t3a(ABCmnK)'};
fun(d3,'D3')

% D4 = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
d4 = {'h2A(ABef),t3a(efCIJK)'};
fun(d4,'D4')

% D5 = einsum_kg(H2A.voov,t3a,'cmke,abeijm->abcijk');
d5 = {'h2A(CmKe),t3a(ABeIJm)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.voov,t3b,'cmke,abeijm->abcijk');
d6 = {'h2B(CmKe),t3b(ABeIJm)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3A,VT3A] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3A_1] = update_t3a_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3A_PPPHHH = %4.15f\n',get_error(X3A(PA,PA,PA,HA,HA,HA),X3A_1))


%% t3A - projection 2 (PPpHHH) [DONE]
clc

out_proj = 'ABcIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
d1 = {'-h1A(mK),t3a(ABcIJm)'};
fun(d1,'D1')

% D2 = einsum_kg(H1A.vv,t3a,'ce,abeijk->abcijk');
d2 = {'h1A(ce),t3a(ABeIJK)'};
fun(d2,'D2')
d3 = {'-h1A(Be),t3a(AceIJK)'};
fun(d3,'D3')

% D3 = 0.5*einsum_kg(H2A.oooo,t3a,'mnij,abcmnk->abcijk');
d4 = {'h2A(mnIJ),t3a(ABcmnK)'};
fun(d4,'D4')

% D4 = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
d5 = {'h2A(ABef),t3a(efcIJK)'};
fun(d5,'D5')
d6 = {'-h2A(Acef),t3a(efBIJK)'};
fun(d6,'D6')

% D5 = einsum_kg(H2A.voov,t3a,'cmke,abeijm->abcijk');
d7 = {'h2A(cmKe),t3a(ABeIJm)'};
fun(d7,'D7')
d8 = {'-h2A(BmKe),t3a(AceIJm)'};
fun(d8,'D8')

% D6 = einsum_kg(H2B.voov,t3b,'cmke,abeijm->abcijk');
d9 = {'h2B(cmKe),t3b(ABeIJm)'};
fun(d9,'D9')
d10 = {'-h2B(BmKe),t3b(AceIJm)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3A,VT3A] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3A_2] = update_t3a_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3A_PPpHHH = %4.15f\n',get_error(X3A(PA,PA,pA,HA,HA,HA),X3A_2))

%% t3A - projection 3 (PPPhHH) [DONE]
clc

out_proj = 'ABCiJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
d1 = {'-h1A(mK),t3a(ABCiJm)'};
fun(d1,'D1')
d2 = {'h1A(mi),t3a(ABCKJm)'};
fun(d2,'D2')

% D2 = einsum_kg(H1A.vv,t3a,'ce,abeijk->abcijk');
d3 = {'h1A(Ce),t3a(ABeiJK)'};
fun(d3,'D3')

% D3 = 0.5*einsum_kg(H2A.oooo,t3a,'mnij,abcmnk->abcijk');
d4 = {'h2A(mniJ),t3a(ABCmnK)'};
fun(d4,'D4')
d5 = {'-h2A(mnKJ),t3a(ABCmni)'};
fun(d5,'D5')

% D4 = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
d6 = {'h2A(ABef),t3a(efCiJK)'};
fun(d6,'D6')

% D5 = einsum_kg(H2A.voov,t3a,'cmke,abeijm->abcijk');
d7 = {'h2A(CmKe),t3a(ABeiJm)'};
fun(d7,'D7')
d8 = {'-h2A(Cmie),t3a(ABeKJm)'};
fun(d8,'D8')

% D6 = einsum_kg(H2B.voov,t3b,'cmke,abeijm->abcijk');
d9 = {'h2B(CmKe),t3b(ABeiJm)'};
fun(d9,'D9')
d10 = {'-h2B(Cmie),t3b(ABeKJm)'};
fun(d10,'D10')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3A,VT3A] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3A_3] = update_t3a_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3A_PPPhHH = %4.15f\n',get_error(X3A(PA,PA,PA,hA,HA,HA),X3A_3))

%% t3A - projection 4 (PPphHH) [DONE]
clc

out_proj = 'ABciJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
d1 = {'-h1A(mK),t3a(ABciJm)'};
fun(d1,'D1')
d2 = {'h1A(mi),t3a(ABcKJm)'};
fun(d2,'D2')

% D2 = einsum_kg(H1A.vv,t3a,'ce,abeijk->abcijk');
d3 = {'h1A(ce),t3a(ABeiJK)'};
fun(d3,'D3')
d4 = {'-h1A(Be),t3a(AceiJK)'};
fun(d4,'D4')

% D3 = 0.5*einsum_kg(H2A.oooo,t3a,'mnij,abcmnk->abcijk');
d5 = {'h2A(mniJ),t3a(ABcmnK)'};
fun(d5,'D5')
d6 = {'-h2A(mnKJ),t3a(ABcmni)'};
fun(d6,'D6')

% D4 = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
d7 = {'h2A(ABef),t3a(efciJK)'};
fun(d7,'D7')
d8 = {'-h2A(Acef),t3a(efBiJK)'};
fun(d8,'D8')

% D5 = einsum_kg(H2A.voov,t3a,'cmke,abeijm->abcijk');
d9 = {'h2A(cmKe),t3a(ABeiJm)'};
fun(d9,'D9')
d10 = {'-h2A(BmKe),t3a(AceiJm)'};
fun(d10,'D10')
d11 = {'-h2A(cmie),t3a(ABeKJm)'};
fun(d11,'D11')
d12 = {'h2A(Bmie),t3a(AceKJm)'};
fun(d12,'D12')

% D6 = einsum_kg(H2B.voov,t3b,'cmke,abeijm->abcijk');
d13 = {'h2B(cmKe),t3b(ABeiJm)'};
fun(d13,'D13')
d14 = {'-h2B(BmKe),t3b(AceiJm)'};
fun(d14,'D14')
d15 = {'-h2B(cmie),t3b(ABeKJm)'};
fun(d15,'D15')
d16 = {'h2B(Bmie),t3b(AceKJm)'};
fun(d16,'D16')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3A,VT3A] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3A_4] = update_t3a_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3A_PPphHH = %4.15f\n',get_error(X3A(PA,PA,pA,hA,HA,HA),X3A_4))

%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end
