clear all
clc
close all

format long

addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests'));

%load emiliano-h2o-1.0-test-ints.mat

load h2o-631g-stretched
%load f2-2.0-pvdz.mat % NEEDS NFZC = 2

%% Set up problem

nfzc = 0; nfzv = 0; 
nact_h = 3; nact_p = 6; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
flag_act_scheme = 0;

sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                           nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);
                       
PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;


                       
%% Full CCSDT

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-11;
ccopts.shift = 0;
flag_full = false;

[cc_t,Ecorr_uccsdt] = uccsdt(sys,ccopts);

% This is important! You need to zero the T3 amps outside of the active
% space before testing them with active space builds. This is because the
% active space builds are coded assuming that all T3 amps outside of the
% space are 0. If they are not 0, the results of diagrams contributing to
% amplitudes within the active space, as computed by complete CCSDT
% updates, can differ!
t3a = zero_t3_outside_act(cc_t.t3a,2,'A',sys); cc_t.t3a = t3a;
t3b = zero_t3_outside_act(cc_t.t3b,2,'B',sys); cc_t.t3b = t3b;
t3c = zero_t3_outside_act(cc_t.t3c,2,'C',sys); cc_t.t3c = t3c;
t3d = zero_t3_outside_act(cc_t.t3d,2,'D',sys); cc_t.t3d = t3d;

t1a = cc_t.t1a; t1b = cc_t.t1b;
t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
t3a = cc_t.t3a; t3b = cc_t.t3b; t3c = cc_t.t3c; t3d = cc_t.t3d;

[HBar_t] = build_ucc_HBar(cc_t,sys,false);
H1A = HBar_t.H1A; H1B = HBar_t.H1B;
H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);

% check active build - test
%[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [t3a_ex,X3A,VT3] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [t3a_1,X3A_1] = update_t3a_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

% [t3a] = zero_t3_outside_act(t3a,2,'A',sys);
% D4_full = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
% D4_full = D4_full - permute(D4_full,[3,2,1,4,5,6]) - permute(D4_full,[1,3,2,4,5,6]);
% 
% D4 = -einsum_kg(H2A.vvvv(PA,PA,PA,pA),T3A.PPpHHH,'ABEf,ECfIJK->ABCIJK')...
% +0.5*einsum_kg(H2A.vvvv(PA,PA,PA,PA),T3A.PPPHHH,'ABEF,EFCIJK->ABCIJK');
% D4 = D4 - permute(D4,[3,2,1,4,5,6]) - permute(D4,[1,3,2,4,5,6]);
% 
% get_error(D4_full(PA,PA,PA,HA,HA,HA),D4)


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

[t3a_ex] = update_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3a_1] = update_t3a_proj1_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3a_PPPHHH = %4.15f\n',get_error(t3a_ex(PA,PA,PA,HA,HA,HA),t3a_1))

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
[t3a_ex] = update_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3a_2] = update_t3a_proj2_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3a_PPpHHH = %4.15f\n',get_error(t3a_ex(PA,PA,pA,HA,HA,HA),t3a_2))

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

t3a_ex = update_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3a_3] = update_t3a_proj3_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3a_PPPhHH = %4.15f\n',get_error(t3a_ex(PA,PA,PA,hA,HA,HA),t3a_3))

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
[t3a_ex] = update_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3a_4] = update_t3a_proj4_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3a_PPPHHH = %4.15f\n',get_error(t3a_ex(PA,PA,pA,hA,HA,HA),t3a_4))

%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end
