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

%% t3B - projection 1 (PPPHHH) [DONE]

clc

out_proj = 'ABCIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
d1 = {'-h1A(MI),t3b(ABCMJK)'};
fun(d1,'D1')

% D2 = -einsum_kg(H1B.oo,t3b,'mk,abcijm->abcijk');
d2 = {'-h1B(mK),t3b(ABCIJm)'};
fun(d2,'D2')

% D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3b(eBCIJK)'};
fun(d3,'D3')

% D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
d4 = {'h1B(Ce),t3b(ABeIJK)'};
fun(d4,'D4')

% D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
d5 = {'h2A(mnIJ),t3b(ABCmnK)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
d6 = {'h2B(mnJK),t3b(ABCImn)'};
fun(d6,'D6')

% D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
d7 = {'h2A(ABef),t3b(efCIJK)'};
fun(d7,'D7')

% D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk'); 
d8 = {'h2B(BCef),t3b(AefIJK)'};
fun(d8,'D8')

% D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');  
d9 = {'h2A(AmIe),t3b(eBCmJK)'};
fun(d9,'D9')

% D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk'); 
d10 = {'h2B(AmIe),t3c(BCeJKm)'};
fun(d10,'D10')

% D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
d11 = {'h2B(mCeK),t3a(ABeIJm)'};
fun(d11,'D11')

% D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
d12 = {'h2C(CmKe),t3b(ABeIJm)'};
fun(d12,'D12')

% D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
d13 = {'-h2B(AmeK),t3b(eBCIJm)'};
fun(d13,'D13')

% D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
d14 = {'-h2B(mCIe),t3b(ABemJK)'};
fun(d14,'D14')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);

[X3B_1] = update_t3b_proj1(cc_t.t1a, cc_t.t1b, cc_t.t2a, cc_t.t2b, cc_t.t2c, T3A, T3B, T3C, T3D, sys, 0.0);

fprintf('\nError in X3B_PPPHHH = %4.15f\n',get_error(X3B(PA,PA,PB,HA,HA,HB),X3B_1))

%% t3B - projection 2 (PPpHHH) [DONE]

clc

out_proj = 'ABcIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3b(ABcmJK)'};
fun(d1,'D1')

% D2 = -einsum_kg(H1B.oo,t3b,'mk,abcijm->abcijk');
d2 = {'-h1B(mK),t3b(ABcIJm)'};
fun(d2,'D2')

% D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3b(eBcIJK)'};
fun(d3,'D3')

% D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
d4 = {'h1B(ce),t3b(ABeIJK)'};
fun(d4,'D4')

% D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
d5 = {'h2A(mnIJ),t3b(ABcmnK)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
d6 = {'h2B(mnJK),t3b(ABcImn)'};
fun(d6,'D6')

% D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
d7 = {'h2A(ABef),t3b(efcIJK)'};
fun(d7,'D7')

% D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk'); 
d8 = {'h2B(Bcef),t3b(AefIJK)'};
fun(d8,'D8')

% D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');  
d9 = {'h2A(AmIe),t3b(eBcmJK)'};
fun(d9,'D9')

% D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk'); 
d10 = {'h2B(AmIe),t3c(BceJKm)'};
fun(d10,'D10')

% D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
d11 = {'h2B(mceK),t3a(ABeIJm)'};
fun(d11,'D11')

% D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
d12 = {'h2C(cmKe),t3b(ABeIJm)'};
fun(d12,'D12')

% D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
d13 = {'-h2B(AmeK),t3b(eBcIJm)'};
fun(d13,'D13')

% D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
d14 = {'-h2B(mcIe),t3b(ABemJK)'};
fun(d14,'D14')


% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_2] = update_t3b_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3B_PPpHHH = %4.15f\n',get_error(X3B(PA,PA,pB,HA,HA,HB),X3B_2))

%% t3B - projection 3 (PpPHHH) [DONE]

clc

out_proj = 'AbCIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3b(AbCmJK)'};
fun(d1,'D1')

% D2 = -einsum_kg(H1B.oo,t3b,'mk,abcijm->abcijk');
d2 = {'-h1B(mK),t3b(AbCIJm)'};
fun(d2,'D2')

% D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3b(ebCIJK)','-h1A(be),t3b(eACIJK)'};
fun(d3,'D3')

% D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
d4 = {'h1B(Ce),t3b(AbeIJK)'};
fun(d4,'D4')

% D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
d5 = {'h2A(mnIJ),t3b(AbCmnK)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
d6 = {'h2B(mnJK),t3b(AbCImn)'};
fun(d6,'D6')

% D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
d7 = {'h2A(Abef),t3b(efCIJK)'};
fun(d7,'D7')

% D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk'); 
d8 = {'h2B(bCef),t3b(AefIJK)','-h2B(ACef),t3b(befIJK)'};
fun(d8,'D8')

% D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');  
d9 = {'h2A(AmIe),t3b(ebCmJK)','-h2A(bmIe),t3b(eACmJK)'};
fun(d9,'D9')

% D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk'); 
d10 = {'h2B(AmIe),t3c(bCeJKm)','-h2B(bmIe),t3c(ACeJKm)'};
fun(d10,'D10')

% D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
d11 = {'h2B(mCeK),t3a(AbeIJm)'};
fun(d11,'D11')

% D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
d12 = {'h2C(CmKe),t3b(AbeIJm)'};
fun(d12,'D12')

% D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
d13 = {'-h2B(AmeK),t3b(ebCIJm)','h2B(bmeK),t3b(eACIJm)'};
fun(d13,'D13')

% D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
d14 = {'-h2B(mCIe),t3b(AbemJK)'};
fun(d14,'D14')


% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_3] = update_t3b_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3B_PpPHHH = %4.15f\n',get_error(X3B(PA,pA,PB,HA,HA,HB),X3B_3))


%% t3B - projection 4 (PPPHHh) [DONE]
clc

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);

[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_4] = update_t3b_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('Error in X3B_PPPHHh = %4.15f\n',get_error(X3B(PA,PA,PB,HA,HA,hB),X3B_4))

%% t3B - projection 5 (PPPhHH) [DONE]
clc

out_proj = 'ABCiJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
d1 = {'-h1A(mi),t3b(ABCmJK)','h1A(mJ),t3b(ABCmiK)'};
fun(d1,'D1')

% D2 = -einsum_kg(H1B.oo,t3b,'mk,abcijm->abcijk');
d2 = {'-h1B(mK),t3b(ABCiJm)'};
fun(d2,'D2')

% D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3b(eBCiJK)'};
fun(d3,'D3')

% D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
d4 = {'h1B(Ce),t3b(ABeiJK)'};
fun(d4,'D4')

% D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
d5 = {'h2A(mniJ),t3b(ABCmnK)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
d6 = {'h2B(mnJK),t3b(ABCimn)','-h2B(mniK),t3b(ABCJmn)'};
fun(d6,'D6')

% D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
d7 = {'h2A(ABef),t3b(efCiJK)'};
fun(d7,'D7')

% D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk'); 
d8 = {'h2B(BCef),t3b(AefiJK)'};
fun(d8,'D8')

% D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');  
d9 = {'h2A(Amie),t3b(eBCmJK)','-h2A(AmJe),t3b(eBCmiK)'};
fun(d9,'D9')

% D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk'); 
d10 = {'h2B(Amie),t3c(BCeJKm)','-h2B(AmJe),t3c(BCeiKm)'};
fun(d10,'D10')

% D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
d11 = {'h2B(mCeK),t3a(ABeiJm)'};
fun(d11,'D11')

% D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
d12 = {'h2C(CmKe),t3b(ABeiJm)'};
fun(d12,'D12')

% D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
d13 = {'-h2B(AmeK),t3b(eBCiJm)'};
fun(d13,'D13')

% D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
d14 = {'-h2B(mCie),t3b(ABemJK)','h2B(mCJe),t3b(ABemiK)'};
fun(d14,'D14')

%
% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_5] = update_t3b_proj5(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3B_PPPhHH = %4.15f\n',get_error(X3B(PA,PA,PB,hA,HA,HB),X3B_5))


%% t3B - projection 6 (PPphHH) [DONE]

clc

out_proj = 'ABciJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);


% h1A(oo)t3b
d1 = {'-h1A(mi),t3b(ABcmJK)','-h1A(mJ),t3b(ABcimK)'};
fun(d1,'D1')

% h1B(oo)t3b
d2 = {'-h1B(mK),t3b(ABciJm)'};
fun(d2,'D2')

% h1A(vv)t3b; A(ab)
% D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3b(eBciJK)'};
fun(d3,'D3')

% h1B(vv)t3b
% D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
d4 = {'h1B(ce),t3b(ABeiJK)'};
fun(d4,'D4')

% h2A(oooo)t3b
% D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
d5 = {'h2A(mniJ),t3b(ABcmnK)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
d6 = {'h2B(mnJK),t3b(ABcimn)','-h2B(mniK),t3b(ABcJmn)'};
fun(d6,'D6')

% D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
d7 = {'h2A(ABef),t3b(efciJK)'};
fun(d7,'D7')

% D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk'); 
d8 = {'h2B(Bcef),t3b(AefiJK)'};
fun(d8,'D8')


% h2A(voov)t3b; A(ab)
d9 = {'h2A(Amie),t3b(eBcmJK)','h2A(AmJe),t3b(eBcimK)'};
fun(d9,'D9')

% h2B(voov)t3c; A(ab)
d10 = {'h2B(Amie),t3c(BceJKm)','-h2B(AmJe),t3c(BceiKm)'};
fun(d10,'D10')

% D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
d11 = {'h2B(mceK),t3a(ABeiJm)'};
fun(d11,'D11')

% D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
d12 = {'h2C(cmKe),t3b(ABeiJm)'};
fun(d12,'D12')

% D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
d13 = {'-h2B(AmeK),t3b(eBciJm)'};
fun(d13,'D13')

% D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
d14 = {'-h2B(mcie),t3b(ABemJK)','h2B(mcJe),t3b(ABemiK)'};
fun(d14,'D14')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_6] = update_t3b_proj6(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3B_PPphHH = %4.15f\n',get_error(X3B(PA,PA,pB,hA,HA,HB),X3B_6))

%% t3B - projection 7 (PPpHHh) [DONE]
clc

out_proj = 'ABcIJk';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3b(ABcmJk)'};
fun(d1,'D1')

% D2 = -einsum_kg(H1B.oo,t3b,'mk,abcijm->abcijk');
d2 = {'-h1B(mk),t3b(ABcIJm)'};
fun(d2,'D2')

% D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3b(eBcIJk)'};
fun(d3,'D3')

% D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
d4 = {'h1B(ce),t3b(ABeIJk)'};
fun(d4,'D4')

% D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
d5 = {'h2A(mnIJ),t3b(ABcmnk)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
d6 = {'h2B(mnJk),t3b(ABcImn)'};
fun(d6,'D6')

% D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
d7 = {'h2A(ABef),t3b(efcIJk)'};
fun(d7,'D7')

% D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk'); 
d8 = {'h2B(Bcef),t3b(AefIJk)'};
fun(d8,'D8')

% D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');  
d9 = {'h2A(AmIe),t3b(eBcmJk)'};
fun(d9,'D9')

% D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk'); 
d10 = {'h2B(AmIe),t3c(BceJkm)'};
fun(d10,'D10')

% D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
d11 = {'h2B(mcek),t3a(ABeIJm)'};
fun(d11,'D11')

% D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
d12 = {'h2C(cmke),t3b(ABeIJm)'};
fun(d12,'D12')

% D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
d13 = {'-h2B(Amek),t3b(eBcIJm)'};
fun(d13,'D13')

% D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
d14 = {'-h2B(mcIe),t3b(ABemJk)'};
fun(d14,'D14')

%
% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_7] = update_t3b_proj7(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3B_PPpHHh = %4.15f\n',get_error(X3B(PA,PA,pB,HA,HA,hB),X3B_7))

%% t3B - projection 8 (PpPhHH) [DONE]
clc

out_proj = 'AbCiJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
d1 = {'-h1A(mi),t3b(AbCmJK)','h1A(mJ),t3b(AbCmiK)'};
fun(d1,'D1')

% D2 = -einsum_kg(H1B.oo,t3b,'mk,abcijm->abcijk');
d2 = {'-h1B(mK),t3b(AbCiJm)'};
fun(d2,'D2')

% D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3b(ebCiJK)','-h1A(be),t3b(eACiJK)'};
fun(d3,'D3')

% D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
d4 = {'h1B(Ce),t3b(AbeiJK)'};
fun(d4,'D4')

% D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
d5 = {'h2A(mniJ),t3b(AbCmnK)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
d6 = {'h2B(mnJK),t3b(AbCimn)','-h2B(mniK),t3b(AbCJmn)'};
fun(d6,'D6')

% D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
d7 = {'h2A(Abef),t3b(efCiJK)'};
fun(d7,'D7')

% D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk'); 
d8 = {'h2B(bCef),t3b(AefiJK)','-h2B(ACef),t3b(befiJK)'};
fun(d8,'D8')

% D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');  
d9 = {'h2A(Amie),t3b(ebCmJK)','-h2A(bmie),t3b(eACmJK)','-h2A(AmJe),t3b(ebCmiK)','h2A(bmJe),t3b(eACmiK)'};
fun(d9,'D9')

% D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk'); 
d10 = {'h2B(Amie),t3c(bCeJKm)','-h2B(bmie),t3c(ACeJKm)','-h2B(AmJe),t3c(bCeiKm)','h2B(bmJe),t3c(ACeiKm)'};
fun(d10,'D10')

% D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
d11 = {'h2B(mCeK),t3a(AbeiJm)'};
fun(d11,'D11')

% D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
d12 = {'h2C(CmKe),t3b(AbeiJm)'};
fun(d12,'D12')

% D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
d13 = {'-h2B(AmeK),t3b(ebCiJm)','h2B(bmeK),t3b(eACiJm)'};
fun(d13,'D13')

% D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
d14 = {'-h2B(mCie),t3b(AbemJK)','h2B(mCJe),t3b(AbemiK)'};
fun(d14,'D14')

%
% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_8] = update_t3b_proj8(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3B_PpPhHH = %4.15f\n',get_error(X3B(PA,pA,PB,hA,HA,HB),X3B_8))


%% t3B - projection 9 (PpPHHh) [DONE]
clc

out_proj = 'AbCIJk';

write_term = @(x,y) write_einsum(x,out_proj,y);

% % diagram 1
% d1 = '-h1a(MI),t3b(AbCMJk)';
% s = ccsdt_act_HBarT3_simple(d1);
% %s = {'-h1a(MI),t3b(AbCMJk)','-h1a(mI),t3b(AbCmJk)'};
% write_term(s,'D1')
% 
% % diagram 2
% d2 = '-h1b(Mk),t3B(AbCIJM)';
% s = ccsdt_act_HBarT3_simple(d2);
% %s = {'-h1b(Mk),t3B(AbCIJM)','-h1b(mk),t3b(AbCIJm)'};
% write_term(s,'D2')
% 
% % diagram 3
% d3 = 'h1a(AE),t3b(EbCIJk)';
% s = ccsdt_act_HBarT3_simple(d3);
% %s = {'h1a(AE),t3b(EbCIJk)','h1a(Ae),t3b(ebCIJk)'};
% write_term(s,'D3')
% 
% % diagram 4
% s = {'h1a(bE),t3b(AECIJk)','h1a(be),t3b(AeCIJk)'};
% write_term(s,'D4')
% 
% % diagram 5
% s = {'h1b(CE),t3b(AbEIJk)','h1B(Ce),t3b(AbeIJk)'};
% write_term(s,'D5')
% 
% % diagram 6
% s = {'h2A(AMIE),t3b(EbCMJk)','h2A(AmIE),t3b(EbCmJk)','h2A(AMIe),t3b(ebCMJk)','h2A(AmIe),t3b(ebCmJk)'};
% write_term(s,'D6')
% 
% % diagram 7
% d7 = 'h2A(bMJE),t3B(AECIMk)';
% s = ccsdt_act_HBarT3_simple(d7);
% %s = {'h2A(bMJE),t3B(AECIMk)','h2a(bmJE),t3b(EACmIk)','h2a(bMJe),t3b(AeCIMk)','-h2a(bmJe),t3b(AeCmIk)'};
% write_term(s,'D7')
% 
% % diagram 8
% s = {'h2b(MCEk),t3a(AEbIMJ)','-h2b(mCEk),t3a(AEbmIJ)','h2b(MCek),t3a(AbeIJM)','-h2b(mCek),t3a(AbemJI)'};
% write_term(s,'D8')
% 
% % diagram 9
% s = {'h2B(AMIE),t3c(bCEJkM)','h2b(AmIE),t3C(bCEJkm)','h2B(AMIe),t3C(bCeJkM)','h2B(AmIe),t3c(bCeJkm)'};
% write_term(s,'D9')
% 
% % diagram 10
% s = {'h2b(bMJE),t3c(ACEIkM)','h2b(bmJE),t3c(ACEIkm)','h2b(bMJe),t3c(ACeIkM)','h2b(bmJe),t3c(ACeIkm)'};
% write_term(s,'D10')
% 
% % diagarm 11
% s = {'h2c(CMkE),t3b(AbEIJM)','h2c(CmkE),t3b(AbEIJm)','h2c(CMke),t3b(AbeIJM)','h2c(Cmke),t3b(AbeIJm)'};
% write_term(s,'D11')
% 
% % diagram 12
% s = {'-h2b(AMEk),t3b(EbCIJM)','-h2B(AmEk),t3b(EbCIJm)','-h2b(AMek),t3b(ebCIJM)','-h2b(Amek),t3b(ebCIJm)'};
% write_term(s,'D12')
% 
% % diagram 13
% s = {'-h2b(bMEk),t3b(AECIJM)','-h2b(bmEk),t3b(AECIJm)','-h2b(bMek),t3b(AeCIJM)','-h2b(bmek),t3b(AeCIJm)'};
% write_term(s,'D13')
% 
% % diagram 14
% s = {'-h2b(MCIE),t3b(AbEMJk)','-h2b(mCIE),t3b(AbEmJk)','-h2b(MCIe),t3b(AbeMJk)','-h2b(mCIe),t3b(AbemJk)'};
% write_term(s,'D14')
% 
% % diagram 15
% s = {'0.5h2a(MNIJ),t3b(AbCMNk)','h2a(mNIJ),t3b(AbCmNk)','0.5h2a(mnIJ),t3b(AbCmnk)'};
% write_term(s,'D15')
% 
% % diagram 16
% s = {'h2b(MNIk),t3b(AbCMJN)','h2b(mNIk),t3b(AbCmJN)','h2b(MnIk),t3b(AbCMJn)','h2b(mnIk),t3b(AbCmJn)'};
% write_term(s,'D16')
% 
% % diagram 17
% s = {'0.5h2a(AbEF),t3b(EFCIJk)','h2a(AbEf),t3b(EfCIJk)','0.5h2a(Abef),t3b(efCIJk)'};
% write_term(s,'D17')
% 
% % diagram 18
% s = {'h2b(ACEF),t3b(EbFIJk)','h2b(ACEf),t3b(EbfIJk)','h2b(ACeF),t3b(ebFIJk)','h2b(ACef),t3b(ebfIJk)'};
% write_term(s,'D18')
% 
% % diagram 19
% s = {'h2b(bCEF),t3b(AEFIJk)','h2b(bCEf),t3b(AEfIJk)','h2b(bCeF),t3b(AeFIJk)','h2b(bCef),t3b(AefIJk)'};
% write_term(s,'D19')

% Vt3B intermediates (vooo)

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_9] = update_t3b_proj9(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3B_PpPHHh = %4.15f\n',get_error(X3B(PA,pA,PB,HA,HA,hB),X3B_9))

%% t3C - projection 1 (PPPHHH) [DONE]
clc

out_proj = 'ABCIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3c(ABCmJK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mJ),t3c(ABCImK)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3c(eBCIJK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(AeCIJK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnJK),t3c(ABCImn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mnIJ),t3c(ABCmnK)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(BCef),t3c(AefIJK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(ABef),t3c(efCIJK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(AmIe),t3c(eBCmJK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(AmIe),t3d(eBCmJK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBeJ),t3b(AeCImK)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(BmJe),t3c(AeCImK)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBIe),t3c(AeCmJK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(AmeJ),t3c(eBCImK)'};
fun(d14,'D14')

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3C_1] = update_t3c_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_PPPHHH = %4.15f\n',get_error(X3C(PA,PB,PB,HA,HB,HB),X3C_1))

%% t3C - projection 2 (PPpHHH) [DONE]
clc

out_proj = 'ABcIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3c(ABcmJK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mJ),t3c(ABcImK)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3c(eBcIJK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(AecIJK)','-h1B(ce),t3c(AeBIJK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnJK),t3c(ABcImn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mnIJ),t3c(ABcmnK)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(Bcef),t3c(AefIJK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(ABef),t3c(efcIJK)','-h2B(Acef),t3c(efBIJK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(AmIe),t3c(eBcmJK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(AmIe),t3d(eBcmJK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBeJ),t3b(AecImK)','-h2B(mceJ),t3b(AeBImK)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(BmJe),t3c(AecImK)','-h2C(cmJe),t3c(AeBImK)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBIe),t3c(AecmJK)','h2B(mcIe),t3c(AeBmJK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(AmeJ),t3c(eBcImK)'};
fun(d14,'D14')

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3C_2] = update_t3c_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_PPpHHH = %4.15f\n',get_error(X3C(PA,PB,pB,HA,HB,HB),X3C_2))

%% t3C - projection 3 (pPPHHH) [DONE]
clc

out_proj = 'aBCIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3c(aBCmJK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mJ),t3c(aBCImK)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(ae),t3c(eBCIJK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(aeCIJK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnJK),t3c(aBCImn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mnIJ),t3c(aBCmnK)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(BCef),t3c(aefIJK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(aBef),t3c(efCIJK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(amIe),t3c(eBCmJK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(amIe),t3d(eBCmJK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBeJ),t3b(aeCImK)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(BmJe),t3c(aeCImK)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBIe),t3c(aeCmJK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(ameJ),t3c(eBCImK)'};
fun(d14,'D14')


[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
%
[X3C_3] = update_t3c_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_pPPHHH = %4.15f\n',get_error(X3C(pA,PB,PB,HA,HB,HB),X3C_3))

%% t3C - projection 4 (PPPHhH) [DONE]
clc

out_proj = 'ABCIjK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3c(ABCmjK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mj),t3c(ABCImK)','h1B(mK),t3c(ABCImj)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3c(eBCIjK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(AeCIjK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnjK),t3c(ABCImn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mnIj),t3c(ABCmnK)','-h2B(mnIK),t3c(ABCmnj)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(BCef),t3c(AefIjK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(ABef),t3c(efCIjK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(AmIe),t3c(eBCmjK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(AmIe),t3d(eBCmjK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBej),t3b(AeCImK)','-h2B(mBeK),t3b(AeCImj)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(Bmje),t3c(AeCImK)','-h2C(BmKe),t3c(AeCImj)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBIe),t3c(AeCmjK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(Amej),t3c(eBCImK)','h2B(AmeK),t3c(eBCImj)'};
fun(d14,'D14')
%
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3C_4] = update_t3c_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_PPPHhH = %4.15f\n',get_error(X3C(PA,PB,PB,HA,hB,HB),X3C_4))


%% t3C - projection 5 (PPPhHH) [DONE]
clc

out_proj = 'ABCiJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mi),t3c(ABCmJK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mJ),t3c(ABCimK)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3c(eBCiJK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(AeCiJK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnJK),t3c(ABCimn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mniJ),t3c(ABCmnK)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(BCef),t3c(AefiJK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(ABef),t3c(efCiJK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(Amie),t3c(eBCmJK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(Amie),t3d(eBCmJK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBeJ),t3b(AeCimK)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(BmJe),t3c(AeCimK)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBie),t3c(AeCmJK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(AmeJ),t3c(eBCimK)'};
fun(d14,'D14')
%
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3C_5] = update_t3c_proj5(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_PPPhHH = %4.15f\n',get_error(X3C(PA,PB,PB,hA,HB,HB),X3C_5))

%% t3C - projection 6 (PPpHhH) [DONE]
clc

out_proj = 'ABcIjK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3c(ABcmjK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mj),t3c(ABcImK)','h1B(mK),t3c(ABcImj)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3c(eBcIjK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(AecIjK)','-h1B(ce),t3c(AeBIjK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnjK),t3c(ABcImn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mnIj),t3c(ABcmnK)','-h2B(mnIK),t3c(ABcmnj)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(Bcef),t3c(AefIjK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(ABef),t3c(efcIjK)','-h2B(Acef),t3c(efBIjK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(AmIe),t3c(eBcmjK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(AmIe),t3d(eBcmjK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBej),t3b(AecImK)','-h2B(mcej),t3b(AeBImK)','-h2B(mBeK),t3b(AecImj)','h2B(mceK),t3b(AeBImj)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(Bmje),t3c(AecImK)','-h2C(cmje),t3c(AeBImK)','-h2C(BmKe),t3c(AecImj)','h2C(cmKe),t3c(AeBImj)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBIe),t3c(AecmjK)','h2B(mcIe),t3c(AeBmjK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(Amej),t3c(eBcImK)','h2B(AmeK),t3c(eBcImj)'};
fun(d14,'D14')

%
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3C_6] = update_t3c_proj6(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_PPpHhH = %4.15f\n',get_error(X3C(PA,PB,pB,HA,hB,HB),X3C_6))

%% t3C - projection 7 (PPphHH) [DONE]
clc

out_proj = 'ABciJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mi),t3c(ABcmJK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mJ),t3c(ABcimK)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3c(eBciJK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(AeciJK)','-h1B(ce),t3c(AeBiJK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnJK),t3c(ABcimn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mniJ),t3c(ABcmnK)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(Bcef),t3c(AefiJK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(ABef),t3c(efciJK)','-h2B(Acef),t3c(efBiJK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(Amie),t3c(eBcmJK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(Amie),t3d(eBcmJK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBeJ),t3b(AecimK)','-h2B(mceJ),t3b(AeBimK)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(BmJe),t3c(AecimK)','-h2C(cmJe),t3c(AeBimK)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBie),t3c(AecmJK)','h2B(mcie),t3c(AeBmJK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(AmeJ),t3c(eBcimK)'};
fun(d14,'D14')

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3C_7] = update_t3c_proj7(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_PPphHH = %4.15f\n',get_error(X3C(PA,PB,pB,hA,HB,HB),X3C_7))

%% t3C - projection 8 (pPPHhH) [DONE]
clc

out_proj = 'aBCIjK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3c(aBCmjK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mj),t3c(aBCImK)','h1B(mK),t3c(aBCImj)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(ae),t3c(eBCIjK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(aeCIjK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnjK),t3c(aBCImn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mnIj),t3c(aBCmnK)','-h2B(mnIK),t3c(aBCmnj)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(BCef),t3c(aefIjK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(aBef),t3c(efCIjK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(amIe),t3c(eBCmjK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(amIe),t3d(eBCmjK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBej),t3b(aeCImK)','-h2B(mBeK),t3b(aeCImj)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(Bmje),t3c(aeCImK)','-h2C(BmKe),t3c(aeCImj)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBIe),t3c(aeCmjK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(amej),t3c(eBCImK)','h2B(ameK),t3c(eBCImj)'};
fun(d14,'D14')

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3C_8] = update_t3c_proj8(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_pPPHhH = %4.15f\n',get_error(X3C(pA,PB,PB,HA,hB,HB),X3C_8))

%% t3C - projection 9 (pPPhHH) [DONE]
clc

out_proj = 'aBCiJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%D1 = -einsum_kg(H1A.oo,t3c,'mi,abcmjk->abcijk');
d1 = {'-h1A(mi),t3c(aBCmJK)'};
fun(d1,'D1')

%D2 = -einsum_kg(H1B.oo,t3c,'mj,abcimk->abcijk');
d2 = {'-h1B(mJ),t3c(aBCimK)'};
fun(d2,'D2')

%D3 = +einsum_kg(H1A.vv,t3c,'ae,ebcijk->abcijk');
d3 = {'h1A(ae),t3c(eBCiJK)'};
fun(d3,'D3')

%D4 = +einsum_kg(H1B.vv,t3c,'be,aecijk->abcijk');
d4 = {'h1B(Be),t3c(aeCiJK)'};
fun(d4,'D4')

%D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');
d5 = {'h2C(mnJK),t3c(aBCimn)'};
fun(d5,'D5')

%D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');
d6 = {'h2B(mniJ),t3c(aBCmnK)'};
fun(d6,'D6')

%D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');
d7 = {'h2C(BCef),t3c(aefiJK)'};
fun(d7,'D7')

%D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');
d8 = {'h2B(aBef),t3c(efCiJK)'};
fun(d8,'D8')

%D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');
d9 = {'h2A(amie),t3c(eBCmJK)'};
fun(d9,'D9')

%D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');
d10 = {'h2B(amie),t3d(eBCmJK)'};
fun(d10,'D10')

%D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');
d11 = {'h2B(mBeJ),t3b(aeCimK)'};
fun(d11,'D11')

%D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');
d12 = {'h2C(BmJe),t3c(aeCimK)'};
fun(d12,'D12')

%D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');
d13 = {'-h2B(mBie),t3c(aeCmJK)'};
fun(d13,'D13')

%D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');
d14 = {'-h2B(ameJ),t3c(eBCimK)'};
fun(d14,'D14')

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3C_9] = update_t3c_proj9(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3C_pPPhHH = %4.15f\n',get_error(X3C(pA,PB,PB,hA,HB,HB),X3C_9))

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

%% Functions + Misc.

% %% Numerical sanity checks for our active space diagrams
% 
% shift = 0.0;
% 
% t1a = cc_t.t1a; t1b = cc_t.t1b;
% t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
% t3a = cc_t.t3a; t3b = cc_t.t3b; t3c = cc_t.t3c; t3d = cc_t.t3d;
% H1A = HBar_t.H1A; H1B = HBar_t.H1B;
% H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
% 
% PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
% PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
% 
% 
% % (Vt3) diagram intermediate - chi(PPvH) (chi(ABeJ))
% D1 = 0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),t3a(PA,PA,PA,HA,HA,HA),'MNeF,ABFMJN->ABeJ');
% D1 = D1 + einsum_kg(sys.vA_oovv(hA,HA,:,PA),t3a(PA,PA,PA,hA,HA,HA),'mNeF,ABFmJN->ABeJ');
% D1 = D1 + 0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),t3a(PA,PA,PA,hA,hA,HA),'mneF,AFBmnJ->ABeJ');
% D1 = D1 + 0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),t3a(PA,PA,pA,HA,HA,HA),'MNef,ABfMJN->ABeJ');
% D1 = D1 + einsum_kg(sys.vA_oovv(hA,HA,:,pA),t3a(PA,PA,pA,hA,HA,HA),'mNef,ABfmJN->ABeJ');
% D1 = D1 - 0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),t3a(PA,PA,pA,hA,hA,HA),'mnef,ABfmnJ->ABeJ');
% 
% D2 = einsum_kg(sys.vB_oovv(HA,HB,:,PB),t3b(PA,PA,PB,HA,HA,HB),'MNeF,ABFMJN->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(hA,HB,:,PB),t3b(PA,PA,PB,hA,HA,HB),'mNeF,ABFmJN->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(HA,hB,:,PB),t3b(PA,PA,PB,HA,HA,hB),'MneF,ABFMJn->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(hA,hB,:,PB),t3b(PA,PA,PB,hA,HA,hB),'mneF,ABFmJn->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(HA,HB,:,pB),t3b(PA,PA,pA,HA,HA,HB),'MNef,ABfMJN->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(hA,HB,:,pB),t3b(PA,PA,pB,hA,HA,HB),'mNef,ABfmJN->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(HA,hB,:,pB),t3b(PA,PA,pB,HA,HA,hB),'Mnef,ABfMJn->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(hA,hA,:,pB),t3b(PA,PA,pB,hA,HA,hB),'mnef,ABfmJn->ABeJ');
% 
% chi2A_PPvH = -D1 - D2;
% 
% X3B_ABCIJk = einsum_kg(chi2A_PPvH,t2b(:,PB,HA,hB),'ABeJ,eCIk->ABCIJk');
% X3B_ABCIJk = X3B_ABCIJk - permute(X3B_ABCIJk,[1,2,3,5,4,6]);
% 
% 
% % 
% chi2A_vvvo = -0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,abfmjn->abej')...
%                  -einsum_kg(sys.vB_oovv,t3b,'mnef,abfmjn->abej');
% 
% X3B_full = einsum_kg(chi2A_vvvo,t2b,'abej,ecik->abcijk');
% X3B_full = X3B_full - permute(X3B_full,[1,2,3,5,4,6]);
% 
% 
% get_error(chi2A_PPvH,chi2A_vvvo(PA,PA,:,HA))
% get_error(X3B_ABCIJk,X3B_full(PA,PA,PB,HA,HA,hB))
% 
% 
% % h2A(oooo) * t3b
% D1 = 0.5*einsum_kg(H2A.oooo(HA,HA,HA,HA),t3b(PA,PA,PB,HA,HA,hB),'MNIJ,ABCMNk->ABCIJk');
% D1 = D1 + einsum_kg(H2A.oooo(hA,HA,HA,HA),t3b(PA,PA,PB,hA,HA,hB),'mNIJ,ABCmNk->ABCIJk');
% D1 = D1 + 0.5*einsum_kg(H2A.oooo(hA,hA,HA,HA),t3b(PA,PA,PB,hA,hA,hB),'mnIJ,ABCmnk->ABCIJk');
% 
% D1_full = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
% 
% get_error(D1,D1_full(PA,PA,PA,HA,HA,hB))
% 
% 
% 
% %% ficticious t vectors
% clear all
% clc
% close all
% 
% load h2o-631g-stretched
% 
% nfzc = 0; nfzv = 0; 
% nact_h = 2; nact_p = 3; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
% nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
% flag_act_scheme = 2;
% 
% sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
%                            nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);
%                        
% PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
% PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;
% 
% t2a = rand(sys.Nvir_alpha,sys.Nvir_alpha,sys.Nocc_alpha,sys.Nocc_alpha);
% t2b = rand(sys.Nvir_alpha,sys.Nvir_beta,sys.Nocc_alpha,sys.Nocc_beta);
% t2c = rand(sys.Nvir_beta,sys.Nvir_beta,sys.Nocc_beta,sys.Nocc_beta);
% 
% t3a = rand(sys.Nvir_alpha,sys.Nvir_alpha,sys.Nvir_alpha,sys.Nocc_alpha,sys.Nocc_alpha,sys.Nocc_alpha);
% t3b = rand(sys.Nvir_alpha,sys.Nvir_alpha,sys.Nvir_beta,sys.Nocc_alpha,sys.Nocc_alpha,sys.Nocc_beta);
% t3c = rand(sys.Nvir_alpha,sys.Nvir_beta,sys.Nvir_beta,sys.Nocc_alpha,sys.Nocc_beta,sys.Nocc_beta);
% t3d = rand(sys.Nvir_beta,sys.Nvir_beta,sys.Nvir_beta,sys.Nocc_beta,sys.Nocc_beta,sys.Nocc_beta);
% 
% t2a = asym(t2a,'a'); t2b = asym(t2b,'b'); t2c = asym(t2c,'c');
% t3a = asym(t3a,'a'); t3b = asym(t3b,'b'); t3c = asym(t3c,'c'); t3d = asym(t3d,'d');
% 
% T3A.PPPHHH = t3a(PA,PA,PA,HA,HA,HA);
% T3A.PPPhHH = t3a(PA,PA,PA,hA,HA,HA);
% T3A.PPPhhH = t3a(PA,PA,PA,hA,hA,HA);
% T3A.PPpHHH = t3a(PA,PA,pA,HA,HA,HA);
% T3A.PPphHH = t3a(PA,PA,pA,hA,HA,HA);
% T3A.PPphhH = t3a(PA,PA,pA,hA,hA,HA);
% 
% T3B.PPPHHH = t3b(PA,PA,PA,HA,HA,HA);
% T3B.PPPhHH = t3b(PA,PA,PA,hA,HA,HA);
% T3B.PPPHHh = t3b(PA,PA,PA,HA,HA,hA);
% T3B.PPPhHh = t3b(PA,PA,PA,hA,HA,hA);
% T3B.PPpHHH = t3b(PA,PA,pA,HA,HA,HA);
% T3B.PPphHH = t3b(PA,PA,pA,hA,HA,HA);
% T3B.PPpHHh = t3b(PA,PA,pA,HA,HA,hA);
% T3B.PPphHh = t3b(PA,PA,pA,hA,HA,hA);
% 
% 
% % (Vt3) diagram intermediate - chi(PPvH) (chi(ABeJ))
% D1 = 0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3A.PPPHHH,'MNeF,ABFMJN->ABeJ');
% D1 = D1 + einsum_kg(sys.vA_oovv(hA,HA,:,PA),T3A.PPPhHH,'mNeF,ABFmJN->ABeJ');
% D1 = D1 + 0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3A.PPPhhH,'mneF,AFBmnJ->ABeJ');
% D1 = D1 + 0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3A.PPpHHH,'MNef,ABfMJN->ABeJ');
% D1 = D1 + einsum_kg(sys.vA_oovv(hA,HA,:,pA),T3A.PPphHH,'mNef,ABfmJN->ABeJ');
% D1 = D1 - 0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3A.PPphhH,'mnef,ABfmnJ->ABeJ');
% 
% D2 = einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3B.PPPHHH,'MNeF,ABFMJN->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3B.PPPhHH,'mNeF,ABFmJN->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3B.PPPHHh,'MneF,ABFMJn->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3B.PPPhHh,'mneF,ABFmJn->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3B.PPpHHH,'MNef,ABfMJN->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3B.PPphHH,'mNef,ABfmJN->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3B.PPpHHh,'Mnef,ABfMJn->ABeJ');
% D2 = D2 + einsum_kg(sys.vB_oovv(hA,hA,:,pB),T3B.PPphHh,'mnef,ABfmJn->ABeJ');
% 
% chi2A_PPvH = -D1 - D2;
% 
% X3B_ABCIJk = einsum_kg(chi2A_PPvH,t2b(:,PB,HA,hB),'ABeJ,eCIk->ABCIJk');
% X3B_ABCIJk = X3B_ABCIJk - permute(X3B_ABCIJk,[1,2,3,5,4,6]);
% 
% 
% % 
% chi2A_vvvo = -0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,abfmjn->abej')...
%                  -einsum_kg(sys.vB_oovv,t3b,'mnef,abfmjn->abej');
% 
% X3B_full = einsum_kg(chi2A_vvvo,t2b,'abej,ecik->abcijk');
% X3B_full = X3B_full - permute(X3B_full,[1,2,3,5,4,6]);
% 
% 
% get_error(chi2A_PPvH,chi2A_vvvo(PA,PA,:,HA))
% get_error(X3B_ABCIJk,X3B_full(PA,PA,PB,HA,HA,hB))
% 
% %%
% clear all
% clc
% close all
% 
% nx = 5;
% ny = 5;
% 
% n = 1000;
% M = n/2;
% 
% lb = -500; ub = 500;
% A = (ub-lb).*rand([nx,n])+lb;
% B = (ub-lb).*rand([n,ny])+lb;
% 
% C = einsum_kg(A,B,'ik,kj->ij');
% 
% C2 =  einsum_kg(A(:,1:M),B(1:M,:),'iK,Kj->ij')...
%      +einsum_kg(A(:,M+1:n),B(M+1:n,:),'ik,kj->ij');
%  
% err = abs(C-C2);
% max(err(:))

function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end

