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

%% t2C - projection 1 (PPHH) [DONE]
clc

out_proj = 'ABIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eABmIJ)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(ABeIJm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(Anef),t3d(eBfIJn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nAfe),t3c(feBnIJ)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnIf),t3d(ABfmJn)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfI),t3c(fABnmJ)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_1] = update_t2c_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_PPHH = %4.15f\n',get_error(X2C(PB,PB,HB,HB),X2C_1(PB,PB,HB,HB)))

%% t2C - projection 2 (PpHH) [DONE]
clc

out_proj = 'AbIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eAbmIJ)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(AbeIJm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(Anef),t3d(ebfIJn)','-h2C(bnef),t3d(eAfIJn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nAfe),t3c(febnIJ)','-h2B(nbfe),t3c(feAnIJ)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnIf),t3d(AbfmJn)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfI),t3c(fAbnmJ)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_2] = update_t2c_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_PpHH = %4.15f\n',get_error(X2C(PB,pB,HB,HB),X2C_2(PB,pB,HB,HB)))


%% t2C - projection 3 (PPhH) [DONE]
clc

out_proj = 'ABiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eABmiJ)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(ABeiJm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(Anef),t3d(eBfiJn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nAfe),t3c(feBniJ)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnif),t3d(ABfmJn)','h2C(mnJf),t3d(ABfmin)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfi),t3c(fABnmJ)','h2B(nmfJ),t3c(fABnmi)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_3] = update_t2c_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_PPhH = %4.15f\n',get_error(X2C(PB,PB,hB,HB),X2C_3(PB,PB,hB,HB)))

%% t2C - projection 4 (PphH) [DONE]
clc

out_proj = 'AbiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eAbmiJ)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(AbeiJm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(Anef),t3d(ebfiJn)','-h2C(bnef),t3d(eAfiJn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nAfe),t3c(febniJ)','-h2B(nbfe),t3c(feAniJ)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnif),t3d(AbfmJn)','h2C(mnJf),t3d(Abfmin)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfi),t3c(fAbnmJ)','h2B(nmfJ),t3c(fAbnmi)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_4] = update_t2c_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_PphH = %4.15f\n',get_error(X2C(PB,pB,hB,HB),X2C_4(PB,pB,hB,HB)))

%% t2C - projection 5 (PPhh) [DONE]
clc

out_proj = 'ABij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eABmij)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(ABeijm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(Anef),t3d(eBfijn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nAfe),t3c(feBnij)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnif),t3d(ABfmjn)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfi),t3c(fABnmj)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_5] = update_t2c_proj5(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_PPhh = %4.15f\n',get_error(X2C(PB,PB,hB,hB),X2C_5(PB,PB,hB,hB)))

%% t2C - projection 6 (ppHH) [DONE]
clc

out_proj = 'abIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eabmIJ)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(abeIJm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(anef),t3d(ebfIJn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nafe),t3c(febnIJ)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnIf),t3d(abfmJn)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfI),t3c(fabnmJ)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_6] = update_t2c_proj6(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_ppHH = %4.15f\n',get_error(X2C(pB,pB,HB,HB),X2C_6(pB,pB,HB,HB)))

%% t2C - projection 7 (Pphh) [DONE]
clc

out_proj = 'Abij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eAbmij)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(Abeijm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(Anef),t3d(ebfijn)','-h2C(bnef),t3d(eAfijn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nAfe),t3c(febnij)','-h2B(nbfe),t3c(feAnij)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnif),t3d(Abfmjn)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfi),t3c(fAbnmj)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_7] = update_t2c_proj7(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_Pphh = %4.15f\n',get_error(X2C(PB,pB,hB,hB),X2C_7(PB,pB,hB,hB)))

%% t2C - projection 8 (pphH) [DONE]
clc

out_proj = 'abiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eabmiJ)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(abeiJm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(anef),t3d(ebfiJn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nafe),t3c(febniJ)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnif),t3d(abfmJn)','h2C(mnJf),t3d(abfmin)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfi),t3c(fabnmJ)','h2B(nmfJ),t3c(fabnmi)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_8] = update_t2c_proj8(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_pphH = %4.15f\n',get_error(X2C(pB,pB,hB,HB),X2C_8(pB,pB,hB,HB)))

%% t2C - projection 9 (pphh) [DONE]
clc

out_proj = 'abij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

%TR_D1 = einsum_kg(h1A_ov,t3c,'me,eabmij->abij');
d1 = {'h1A(me),t3c(eabmij)'};
fun(d1,'D1')

%TR_D2 = einsum_kg(h1B_ov,t3d,'me,abeijm->abij');
d2 = {'h1B(me),t3d(abeijm)'};
fun(d2,'D2')

%TR_D3 = 0.5*einsum_kg(h2C_vovv,t3d,'anef,ebfijn->abij');
d3 = {'h2C(anef),t3d(ebfijn)'};
fun(d3,'D3')

%TR_D4 = einsum_kg(h2B_ovvv,t3c,'nafe,febnij->abij');
d4 = {'h2B(nafe),t3c(febnij)'};
fun(d4,'D4')

%TR_D5 = -0.5*einsum_kg(h2C_ooov,t3d,'mnif,abfmjn->abij');
d5 = {'-h2C(mnif),t3d(abfmjn)'};
fun(d5,'D5')

%TR_D6 = -einsum_kg(h2B_oovo,t3c,'nmfi,fabnmj->abij');
d6 = {'-h2B(nmfi),t3c(fabnmj)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2C] = build_t2c(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2C_9] = update_t2c_proj9(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2C_pphh = %4.15f\n',get_error(X2C(pB,pB,hB,hB),X2C_9(pB,pB,hB,hB)))

%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end

