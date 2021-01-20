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

% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_1] = update_t3c_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_PPPHHH = %4.15f\n',get_error(X3C(PA,PB,PB,HA,HB,HB),X3C_1))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_1] = update_t3c_proj1_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_PPPHHH = %4.15f\n',get_error(t3c_ex(PA,PB,PB,HA,HB,HB),t3c_1))

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

% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_2] = update_t3c_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_PPpHHH = %4.15f\n',get_error(X3C(PA,PB,pB,HA,HB,HB),X3C_2))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_2] = update_t3c_proj2_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_PPpHHH = %4.15f\n',get_error(t3c_ex(PA,PB,pB,HA,HB,HB),t3c_2))

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


% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_3] = update_t3c_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_pPPHHH = %4.15f\n',get_error(X3C(pA,PB,PB,HA,HB,HB),X3C_3))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_3] = update_t3c_proj3_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_pPPHHH = %4.15f\n',get_error(t3c_ex(pA,PB,PB,HA,HB,HB),t3c_3))

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
% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_4] = update_t3c_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_PPPHhH = %4.15f\n',get_error(X3C(PA,PB,PB,HA,hB,HB),X3C_4))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_4] = update_t3c_proj4_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_PPPHhH = %4.15f\n',get_error(t3c_ex(PA,PB,PB,HA,hB,HB),t3c_4))

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
% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_5] = update_t3c_proj5(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_PPPhHH = %4.15f\n',get_error(X3C(PA,PB,PB,hA,HB,HB),X3C_5))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_5] = update_t3c_proj5_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_PPPhHH = %4.15f\n',get_error(t3c_ex(PA,PB,PB,hA,HB,HB),t3c_5))

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
% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_6] = update_t3c_proj6(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_PPpHhH = %4.15f\n',get_error(X3C(PA,PB,pB,HA,hB,HB),X3C_6))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_6] = update_t3c_proj6_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_PPpHhH = %4.15f\n',get_error(t3c_ex(PA,PB,pB,HA,hB,HB),t3c_6))

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

% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_7] = update_t3c_proj7(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_PPphHH = %4.15f\n',get_error(X3C(PA,PB,pB,hA,HB,HB),X3C_7))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_7] = update_t3c_proj7_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_PPphHH = %4.15f\n',get_error(t3c_ex(PA,PB,pB,hA,HB,HB),t3c_7))

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

% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_8] = update_t3c_proj8(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_pPPHhH = %4.15f\n',get_error(X3C(pA,PB,PB,HA,hB,HB),X3C_8))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_8] = update_t3c_proj8_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_pPPHhH = %4.15f\n',get_error(t3c_ex(pA,PB,PB,HA,hB,HB),t3c_8))

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

% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% [X3C,VT3B,VT3C] = build_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
% [X3C_9] = update_t3c_proj9(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
% fprintf('\nError in X3C_pPPhHH = %4.15f\n',get_error(X3C(pA,PB,PB,hA,HB,HB),X3C_9))

[t3c_ex] = update_t3c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3c_9] = update_t3c_proj9_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,HBar_t,sys,0.0);
fprintf('\nError in t3c_pPPhHH = %4.15f\n',get_error(t3c_ex(pA,PB,PB,hA,HB,HB),t3c_9))

%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end
