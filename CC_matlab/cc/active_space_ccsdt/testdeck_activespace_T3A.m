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
nact_h = 2; nact_p = 4; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
flag_act_scheme = 0;

sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                           nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);
                       
PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;


                       
%% Full CCSDT

ccopts.diis_size = 5;
ccopts.maxit = 3;
ccopts.tol = 1e-11;
ccopts.shift = 0;
flag_full = false;

[cc_t,Ecorr_uccsdt] = uccsdt(sys,ccopts);

% t1a = 100*cc_t.t1a; t1b = 100*cc_t.t1b;
% t2a = 1000*cc_t.t2a; t2b = 1000*cc_t.t2b; t2c = 1000*cc_t.t2c;
% t3a = 10000*cc_t.t3a; t3b = 10000*cc_t.t3b; t3c = 10000*cc_t.t3c; t3d = 10000*cc_t.t3d;

t1a = cc_t.t1a; t1b = cc_t.t1b;
t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
t3a = cc_t.t3a; t3b = cc_t.t3b; t3c = cc_t.t3c; t3d = cc_t.t3d;

% cc_t.t1a = t1a; cc_t.t1b = t1b;
% cc_t.t2a = t2a; cc_t.t2b = t2b; cc_t.t2c = t2c;
% cc_t.t3a = t3a; cc_t.t3b = t3b; cc_t.t3c = t3c; cc_t.t3d = t3d;

[HBar_t] = build_ucc_HBar(cc_t,sys,false);

H1A = HBar_t.H1A; H1B = HBar_t.H1B;
H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

% %% Ficticious T vectors
% 
% Noa  = sys.Nocc_alpha; Nua = sys.Nvir_alpha;
% Nob = sys.Nocc_beta;   Nub = sys.Nvir_beta;
% 
% t1a = rand(Nua,Noa);
% t1b = rand(Nub,Nob);
% t2a = rand(Nua,Nua,Noa,Noa);
% t2b = rand(Nua,Nub,Noa,Nob);
% t2c = rand(Nub,Nub,Nob,Nob);
% t3a = rand(Nua,Nua,Nua,Noa,Noa,Noa);
% t3b = rand(Nua,Nua,Nub,Noa,Noa,Nob);
% t3c = rand(Nua,Nub,Nub,Noa,Nob,Nob);
% t3d = rand(Nub,Nub,Nub,Nob,Nob,Nob);
% 
% t2a = 1/4*(t2a - permute(t2a,[2,1,3,4]) - permute(t2a,[1,2,4,3]) + permute(t2a,[2,1,4,3]));
% t2c = 1/4*(t2c - permute(t2c,[2,1,3,4]) - permute(t2c,[1,2,4,3]) + permute(t2c,[2,1,4,3]));
% 
% t3a = 1/36*(t3a - permute(t3a,[2,1,3,4,5,6]) - permute(t3a,[1,3,2,4,5,6]) ...
%                 - permute(t3a,[3,2,1,4,5,6]) + permute(t3a,[2,3,1,4,5,6]) + permute(t3a,[3,1,2,4,5,6]) ...
%                 - permute(t3a,[1,2,3,5,4,6]) - permute(t3a,[1,2,3,4,6,5]) - permute(t3a,[1,2,3,
% 
% cc_t.t1a = t1a; cc_t.t1b = t1b;
% cc_t.t2a = t2a; cc_t.t2b = t2b; cc_t.t2c = t2c;
% cc_t.t3a = t3a; cc_t.t3b = t3b; cc_t.t3c = t3c; cc_t.t3d = t3d;
% 
% [HBar_t] = build_ucc_HBar(cc_t,sys,false);
% 
% H1A = HBar_t.H1A; H1B = HBar_t.H1B;
% H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
% %%
% clc
% 
% [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
% 
% % active space scheme
% D1 = +einsum_kg(H1A.oo(hA,HA),T3A.PPphhH,'mK,ABcimJ->ABciJK')...
% -einsum_kg(H1A.oo(HA,HA),T3A.PPphHH,'MK,ABciJM->ABciJK');
% D1 = D1 - permute(D1,[1,2,3,4,6,5]);
% 
% D2 = -einsum_kg(H1A.oo(hA,hA),T3A.PPphHH,'mi,ABcmJK->ABciJK')...
% +einsum_kg(H1A.oo(HA,hA),T3A.PPpHHH,'Mi,ABcKJM->ABciJK');
% D1 = D1 + D2;
% 
% % full T3
% D1_full = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
% D1_full = D1_full - permute(D1_full,[1,2,3,4,6,5]) - permute(D1_full,[1,2,3,6,5,4]);
% D1_sub = D1_full(PA,PA,pA,hA,HA,HA);
% 
% get_error(D1,D1_sub)
% 
% %%
% D1_2 = zeros(size(D1));
% for a = ia
%     for b = ib
%         for c = ic
%             for i = ii
%                 for j = ij
%                     for k = ik
%                         tmp1 = 0.0; tmp2 = 0.0;
%                         for m = HA
%                             val1 = -H1A.oo(m,k)*t3a(a,b,c,i,j,m);
%                             val2 = H1A.oo(m,j)*t3a(a,b,c,i,k,m);
%                             val3 = H1A.oo(m,i)*t3a(a,b,c,k,j,m);
%                             tmp1 = tmp1 + val1 + val2 + val3;
%                         end
%                         for m = hA
%                             val1 = -H1A.oo(m,k)*t3a(a,b,c,i,j,m);
%                             val2 = H1A.oo(m,j)*t3a(a,b,c,i,k,m);
%                             val3 = H1A.oo(m,i)*t3a(a,b,c,k,j,m);
%                             tmp2 = tmp2 + val1 + val2 + val3;
%                         end
%                         D1_2(a,b,c,i,j,k) = tmp1 + tmp2;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% tmp1 + tmp2
% 
% % unravel full T3 in loops
% D1_full2 = zeros(size(D1_full));
% for a = ia
%     for b = ib
%         for c = ic
%             for i = ii
%                 for j = ij
%                     for k = ik
%                         tmp1 = 0.0; tmp2 = 0.0;
%                         for m = HA
%                             val1 = -H1A.oo(m,k)*t3a(a,b,c,i,j,m);
%                             val2 = H1A.oo(m,j)*t3a(a,b,c,i,k,m);
%                             val3 = H1A.oo(m,i)*t3a(a,b,c,k,j,m);
%                             tmp1 = tmp1 + val1 + val2 + val3;
%                         end
%                         for m = hA
%                             val1 = -H1A.oo(m,k)*t3a(a,b,c,i,j,m);
%                             val2 = H1A.oo(m,j)*t3a(a,b,c,i,k,m);
%                             val3 = H1A.oo(m,i)*t3a(a,b,c,k,j,m);
%                             tmp2 = tmp2 + val1 + val2 + val3;
%                         end
%                         D1_full2(a,b,c,i,j,k) = tmp1 + tmp2;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% tmp1 + tmp2
% %%
% D1_sub2 = D1_full2(PA,PA,pA,hA,HA,HA);
% 
% get_error(D1_sub,D1_sub2)
% get_error(D1,D1_sub2)

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
%[t3a_ex,X3A,VT3] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
%[t3a_1,X3A_1] = update_t3a_proj1(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0)

t3a_ex = update_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3a_ex] = zero_t3_outside_act(t3a_ex,2,'A',sys);
[t3a_1] = update_t3a_proj1_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

%fprintf('\nError in X3A_PPPHHH = %4.15f\n',get_error(X3A(PA,PA,PA,HA,HA,HA),X3A_1))
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
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[t3a_ex,X3A,VT3A] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3a_2,X3A_2] = update_t3a_proj2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3A_PPpHHH = %4.15f\n',get_error(X3A(PA,PA,pA,HA,HA,HA),X3A_2))
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

%%
clc

%[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);

cc_t.t1a = t1a; cc_t.t1b = t1b; cc_t.t2a = t2a; cc_t.t2b = t2b; cc_t.t2c =  t2c;
[HBar_t] = build_ucc_HBar(cc_t,sys,false);

H1A = HBar_t.H1A; H1B = HBar_t.H1B;
H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

D1 = +einsum_kg(H1A.oo(hA,HA),T3A.PPPhhH,'mK,ABCimJ->ABCiJK')...
-einsum_kg(H1A.oo(HA,HA),T3A.PPPhHH,'MK,ABCiJM->ABCiJK');
D1 = D1 - permute(D1,[1,2,3,4,6,5]);

D2 = +einsum_kg(H1A.oo(hA,hA),T3A.PPPhHH,'mi,ABCmKJ->ABCiJK')...
+einsum_kg(H1A.oo(HA,hA),T3A.PPPHHH,'Mi,ABCKJM->ABCiJK');

D3 = +einsum_kg(H1A.vv(PA,pA),T3A.PPphHH,'Ce,ABeiJK->ABCiJK')...
+einsum_kg(H1A.vv(PA,PA),T3A.PPPhHH,'CE,ABEiJK->ABCiJK');
D3 = D3 - permute(D3,[3,2,1,4,5,6]) - permute(D3,[1,3,2,4,5,6]);

D4 = +0.5*einsum_kg(H2A.oooo(hA,hA,hA,HA),T3A.PPPhhH,'mniJ,ABCmnK->ABCiJK')...
-einsum_kg(H2A.oooo(HA,hA,hA,HA),T3A.PPPhHH,'MniJ,ABCnMK->ABCiJK')...
+0.5*einsum_kg(H2A.oooo(HA,HA,hA,HA),T3A.PPPHHH,'MNiJ,ABCMNK->ABCiJK');
D4 = D4 - permute(D4,[1,2,3,4,6,5]);

D5 = -0.5*einsum_kg(H2A.oooo(hA,hA,HA,HA),T3A.PPPhhh,'mnKJ,ABCmni->ABCiJK')...
-einsum_kg(H2A.oooo(HA,hA,HA,HA),T3A.PPPhhH,'MnKJ,ABCniM->ABCiJK')...
-0.5*einsum_kg(H2A.oooo(HA,HA,HA,HA),T3A.PPPhHH,'MNKJ,ABCiMN->ABCiJK');

D6 = +0.5*einsum_kg(H2A.vvvv(PA,PA,pA,pA),T3A.PpphHH,'ABef,CefiJK->ABCiJK')...
-einsum_kg(H2A.vvvv(PA,PA,PA,pA),T3A.PPphHH,'ABEf,ECfiJK->ABCiJK')...
+0.5*einsum_kg(H2A.vvvv(PA,PA,PA,PA),T3A.PPPhHH,'ABEF,EFCiJK->ABCiJK');
D6 = D6 - permute(D6,[3,2,1,4,5,6]) - permute(D6,[1,3,2,4,5,6]);

D7 = -einsum_kg(H2A.voov(PA,hA,HA,pA),T3A.PPphhH,'CmKe,ABeimJ->ABCiJK')...
-einsum_kg(H2A.voov(PA,hA,HA,PA),T3A.PPPhhH,'CmKE,ABEimJ->ABCiJK')...
+einsum_kg(H2A.voov(PA,HA,HA,pA),T3A.PPphHH,'CMKe,ABeiJM->ABCiJK')...
+einsum_kg(H2A.voov(PA,HA,HA,PA),T3A.PPPhHH,'CMKE,ABEiJM->ABCiJK');
D7 = D7 - permute(D7,[1,2,3,4,6,5]);

D8 = -einsum_kg(H2A.voov(PA,hA,hA,pA),T3A.PPphHH,'Cmie,ABemKJ->ABCiJK')...
-einsum_kg(H2A.voov(PA,hA,hA,PA),T3A.PPPhHH,'CmiE,ABEmKJ->ABCiJK')...
-einsum_kg(H2A.voov(PA,HA,hA,pA),T3A.PPpHHH,'CMie,ABeKJM->ABCiJK')...
-einsum_kg(H2A.voov(PA,HA,hA,PA),T3A.PPPHHH,'CMiE,ABEKJM->ABCiJK');

D78 = D7 + D8;
D78 = D78 - permute(D78,[3,2,1,4,5,6]) - permute(D78,[1,3,2,4,5,6]);

D9 = +einsum_kg(H2B.voov(PA,hB,HA,pB),T3B.PPphHh,'CmKe,ABeiJm->ABCiJK')...
+einsum_kg(H2B.voov(PA,hB,HA,PB),T3B.PPPhHh,'CmKE,ABEiJm->ABCiJK')...
+einsum_kg(H2B.voov(PA,HB,HA,pB),T3B.PPphHH,'CMKe,ABeiJM->ABCiJK')...
+einsum_kg(H2B.voov(PA,HB,HA,PB),T3B.PPPhHH,'CMKE,ABEiJM->ABCiJK');
D9 = D9 - permute(D9,[1,2,3,4,6,5]);

D10 = -einsum_kg(H2B.voov(PA,hB,hA,pB),T3B.PPpHHh,'Cmie,ABeKJm->ABCiJK')...
-einsum_kg(H2B.voov(PA,hB,hA,PB),T3B.PPPHHh,'CmiE,ABEKJm->ABCiJK')...
-einsum_kg(H2B.voov(PA,HB,hA,pB),T3B.PPpHHH,'CMie,ABeKJM->ABCiJK')...
-einsum_kg(H2B.voov(PA,HB,hA,PB),T3B.PPPHHH,'CMiE,ABEKJM->ABCiJK');

D910 = D9 + D10;
D910 = D910 - permute(D910,[3,2,1,4,5,6]) - permute(D910,[1,3,2,4,5,6]);

%
D1_full = -einsum_kg(H1A.oo,t3a,'mk,abcijm->abcijk');
D1_full = D1_full - permute(D1_full,[1,2,3,6,5,4]) - permute(D1_full,[1,2,3,4,6,5]);
D2_full = einsum_kg(H1A.vv,t3a,'ce,abeijk->abcijk');
D2_full = D2_full - permute(D2_full,[3,2,1,4,5,6]) - permute(D2_full,[1,3,2,4,5,6]);
D3_full = 0.5*einsum_kg(H2A.oooo,t3a,'mnij,abcmnk->abcijk'); 
D3_full = D3_full - permute(D3_full,[1,2,3,6,5,4]) - permute(D3_full,[1,2,3,4,6,5]);
D4_full = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
D4_full = D4_full - permute(D4_full,[3,2,1,4,5,6]) - permute(D4_full,[1,3,2,4,5,6]);
D5_full = einsum_kg(H2A.voov,t3a,'cmke,abeijm->abcijk');
D5_full = D5_full - permute(D5_full,[1,2,3,6,5,4]) - permute(D5_full,[1,2,3,4,6,5])...
        - permute(D5_full,[3,2,1,4,5,6]) - permute(D5_full,[1,3,2,4,5,6]) ...
        + permute(D5_full,[3,2,1,6,5,4]) + permute(D5_full,[3,2,1,4,6,5]) ...
        + permute(D5_full,[1,3,2,6,5,4]) + permute(D5_full,[1,3,2,4,6,5]); 
D6_full = einsum_kg(H2B.voov,t3b,'cmke,abeijm->abcijk');
D6_full = D6_full - permute(D6_full,[1,2,3,6,5,4]) - permute(D6_full,[1,2,3,4,6,5])...
        - permute(D6_full,[3,2,1,4,5,6]) - permute(D6_full,[1,3,2,4,5,6]) ...
        + permute(D6_full,[3,2,1,6,5,4]) + permute(D6_full,[3,2,1,4,6,5]) ...
        + permute(D6_full,[1,3,2,6,5,4]) + permute(D6_full,[1,3,2,4,6,5]); 


X3A_3 = D1 + D2 + D3 + D4 + D5 + D6 + D78 + D910;
X3A = D1_full + D2_full;% + D3_full + D4_full + D5_full + D6_full;
X3A = zero_t3_outside_act(X3A,2,'A',sys);

fprintf('\nError in X3A_PPPhHH = %4.15f\n',get_error(X3A(PA,PA,PA,hA,HA,HA),X3A_3))
%%
clc

% check active build
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[t3a_ex,X3A,VT3A] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3a_3,X3A_3] = update_t3a_proj3(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);
%t3a_ex = update_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
%[t3a_ex] = zero_t3_outside_act(t3a_ex,2,'A',sys);
%[t3a_3] = update_t3a_proj3_ccsdt2_v2(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3A_PPPhHH = %4.15f\n',get_error(X3A(PA,PA,PA,hA,HA,HA),X3A_3))
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
[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
[t3a_ex,X3A,VT3A] = build_t3a(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[t3a_4,X3A_4] = update_t3a_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('\nError in X3A_PPphHH = %4.15f\n',get_error(X3A(PA,PA,pA,hA,HA,HA),X3A_4))
fprintf('\nError in t3a_PPphHH = %4.15f\n',get_error(t3a_ex(PA,PA,pA,hA,HA,HA),t3a_4))

%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end
