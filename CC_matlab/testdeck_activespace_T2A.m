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

%% t2A - projection 1 (PPHH) [DONE]
clc

out_proj = 'ABIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(ABeIJm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(ABeIJm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnIf),t3b(ABfmJn)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnIf),t3a(ABfmJn)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(Anef),t3a(eBfIJn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(Anef),t3b(eBfIJn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_1] = update_t2a_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_PPHH = %4.15f\n',get_error(X2A(PA,PA,HA,HA),X2A_1(PA,PA,HA,HA)))

%% t2A - projection 2 (PpHH) [DONE]
clc

out_proj = 'AbIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(AbeIJm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(AbeIJm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnIf),t3b(AbfmJn)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnIf),t3a(AbfmJn)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(Anef),t3a(ebfIJn)','-h2A(bnef),t3a(eAfIJn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(Anef),t3b(ebfIJn)','-h2B(bnef),t3b(eAfIJn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_2] = update_t2a_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_PpHH = %4.15f\n',get_error(X2A(PA,pA,HA,HA),X2A_2(PA,pA,HA,HA)))

%% t2A - projection 3 (PPhH) 
clc

out_proj = 'ABiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(ABeiJm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(ABeiJm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnif),t3b(ABfmJn)','h2B(mnJf),t3b(ABfmin)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnif),t3a(ABfmJn)','h2A(mnJf),t3a(ABfmin)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(Anef),t3a(eBfiJn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(Anef),t3b(eBfiJn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_3] = update_t2a_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_PPhH = %4.15f\n',get_error(X2A(PA,PA,hA,HA),X2A_3(PA,PA,hA,HA)))

%% t2A - projection 4 (PphH) [DONE]
clc

out_proj = 'AbiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(AbeiJm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(AbeiJm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnif),t3b(AbfmJn)','h2B(mnJf),t3b(Abfmin)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnif),t3a(AbfmJn)','h2A(mnJf),t3a(Abfmin)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(Anef),t3a(ebfiJn)','-h2A(bnef),t3a(eAfiJn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(Anef),t3b(ebfiJn)','-h2B(bnef),t3b(eAfiJn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_4] = update_t2a_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_PphH = %4.15f\n',get_error(X2A(PA,pA,hA,HA),X2A_4(PA,pA,hA,HA)))

%% t2A - projection 5 (PPhh) [DONE]
clc

out_proj = 'ABij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(ABeijm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(ABeijm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnif),t3b(ABfmjn)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnif),t3a(ABfmjn)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(Anef),t3a(eBfijn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(Anef),t3b(eBfijn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_5] = update_t2a_proj5(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_PPhh = %4.15f\n',get_error(X2A(PA,PA,hA,hA),X2A_5(PA,PA,hA,hA)))

%% t2A - projection 6 (ppHH) [DONE]
clc

out_proj = 'abIJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(abeIJm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(abeIJm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnIf),t3b(abfmJn)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnIf),t3a(abfmJn)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(anef),t3a(ebfIJn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(anef),t3b(ebfIJn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_6] = update_t2a_proj6(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_ppHH = %4.15f\n',get_error(X2A(pA,pA,HA,HA),X2A_6(pA,pA,HA,HA)))

%% t2A - projection 7 (Pphh) [DONE]
clc

out_proj = 'Abij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(Abeijm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(Abeijm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnif),t3b(Abfmjn)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnif),t3a(Abfmjn)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(Anef),t3a(ebfijn)','-h2A(bnef),t3a(eAfijn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(Anef),t3b(ebfijn)','-h2B(bnef),t3b(eAfijn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_7] = update_t2a_proj7(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_Pphh = %4.15f\n',get_error(X2A(PA,pA,hA,hA),X2A_7(PA,pA,hA,hA)))

%% t2A - projection 8 (pphH) [DONE]
clc

out_proj = 'abiJ';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(abeiJm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(abeiJm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnif),t3b(abfmJn)','h2B(mnJf),t3b(abfmin)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnif),t3a(abfmJn)','h2A(mnJf),t3a(abfmin)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(anef),t3a(ebfiJn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(anef),t3b(ebfiJn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_8] = update_t2a_proj8(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_pphH = %4.15f\n',get_error(X2A(pA,pA,hA,HA),X2A_8(pA,pA,hA,HA)))

%% t2A - projection 9 (pphh) [DONE]
clc

out_proj = 'abij';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% TR_D1 = einsum_kg(h1A_ov,t3a,'me,abeijm->abij');
d1 = {'h1A(me),t3a(abeijm)'};
fun(d1,'D1')

% TR_D2 = einsum_kg(h1B_ov,t3b,'me,abeijm->abij');
d2 = {'h1B(me),t3b(abeijm)'};
fun(d2,'D2')
 
% TR_D3 = -einsum_kg(h2B_ooov,t3b,'mnif,abfmjn->abij');
d3 = {'-h2B(mnif),t3b(abfmjn)'};
fun(d3,'D3')
 
% TR_D4 = -0.5*einsum_kg(h2A_ooov,t3a,'mnif,abfmjn->abij');
d4 = {'-h2A(mnif),t3a(abfmjn)'};
fun(d4,'D4')

% TR_D5 = 0.5*einsum_kg(h2A_vovv,t3a,'anef,ebfijn->abij');
d5 = {'h2A(anef),t3a(ebfijn)'};
fun(d5,'D5')
 
% TR_D6 = einsum_kg(h2B_vovv,t3b,'anef,ebfijn->abij');
d6 = {'h2B(anef),t3b(ebfijn)'};
fun(d6,'D6')

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[X2A] = build_t2a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X2A_9] = update_t2a_proj9(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X2A_pphh = %4.15f\n',get_error(X2A(pA,pA,hA,hA),X2A_9(pA,pA,hA,hA)))
%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end

