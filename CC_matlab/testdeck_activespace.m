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
nact_h = 2; nact_p = 3; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
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

% HBar
flag_build_3body = true;
[HBar_t] = build_ucc_HBar( cc_t, sys, flag_build_3body);

%%

t1a = cc_t.t1a; t1b = cc_t.t1b;
t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
t3a = cc_t.t3a; t3b = cc_t.t3b; t3c = cc_t.t3c; t3d = cc_t.t3d;
H1A = HBar_t.H1A; H1B = HBar_t.H1B;
H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3C = HBar_t.H3C; H3D = HBar_t.H3D;


A1 = einsum_kg(H2A.ooov,t2a,'mnie,abnj->abmije');
A1 = A1 - permute(A1,[1,2,3,5,4,6]);

A2 = einsum_kg(H2A.vovv,t2a,'amfe,fbij->abmije');
A2 = A2 - permute(A2,[2,1,3,4,5,6]);

X = A1 + A2;

get_error(X,H3A.vvooov)

%% t3B - projection 4 (PPPHHh)
clc

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);

[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
[X3B_4,Vt3A,Vt3B] = update_t3b_proj4(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,T3A,T3B,T3C,T3D,sys,0.0);

fprintf('Error in X3B_PPPHHh = %4.15f\n',get_error(X3B(PA,PA,PB,HA,HA,hB),X3B_4))
%fprintf('Error in Vt3A_oPHH = %4.15f\n',get_error(VT3A.vooo(PA,:,HA,HA),Vt3A.PoHH))
%fprintf('Error in Vt3B_PoHH = %4.15f\n',get_error(VT3B.vooo(PA,:,HA,hB),Vt3B.PoHh))
%fprintf('Error in Vt3B_oPHh = %4.15f\n',get_error(VT3B.ovoo(:,PB,HA,hB),Vt3B.oPHh))
%fprintf('Error in Vt3A_PPHv = %4.15f\n',get_error(VT3A.vvov(PA,PA,HA,:),Vt3A.PPHv))
%fprintf('Error in Vt3B_PPvh = %4.15f\n',get_error(VT3B.vvvo(PA,PB,:,hB),Vt3B.PPvh))
%fprintf('Error in Vt3B_PPHv = %4.15f\n',get_error(VT3B.vvov(PA,PB,HA,:),Vt3B.PPHv))

%% Numerical sanity checks for our active space diagrams

shift = 0.0;

t1a = cc_t.t1a; t1b = cc_t.t1b;
t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
t3a = cc_t.t3a; t3b = cc_t.t3b; t3c = cc_t.t3c; t3d = cc_t.t3d;
H1A = HBar_t.H1A; H1B = HBar_t.H1B;
H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;


% (Vt3) diagram intermediate - chi(PPvH) (chi(ABeJ))
D1 = 0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),t3a(PA,PA,PA,HA,HA,HA),'MNeF,ABFMJN->ABeJ');
D1 = D1 + einsum_kg(sys.vA_oovv(hA,HA,:,PA),t3a(PA,PA,PA,hA,HA,HA),'mNeF,ABFmJN->ABeJ');
D1 = D1 + 0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),t3a(PA,PA,PA,hA,hA,HA),'mneF,AFBmnJ->ABeJ');
D1 = D1 + 0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),t3a(PA,PA,pA,HA,HA,HA),'MNef,ABfMJN->ABeJ');
D1 = D1 + einsum_kg(sys.vA_oovv(hA,HA,:,pA),t3a(PA,PA,pA,hA,HA,HA),'mNef,ABfmJN->ABeJ');
D1 = D1 - 0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),t3a(PA,PA,pA,hA,hA,HA),'mnef,ABfmnJ->ABeJ');

D2 = einsum_kg(sys.vB_oovv(HA,HB,:,PB),t3b(PA,PA,PB,HA,HA,HB),'MNeF,ABFMJN->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(hA,HB,:,PB),t3b(PA,PA,PB,hA,HA,HB),'mNeF,ABFmJN->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(HA,hB,:,PB),t3b(PA,PA,PB,HA,HA,hB),'MneF,ABFMJn->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(hA,hB,:,PB),t3b(PA,PA,PB,hA,HA,hB),'mneF,ABFmJn->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(HA,HB,:,pB),t3b(PA,PA,pA,HA,HA,HB),'MNef,ABfMJN->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(hA,HB,:,pB),t3b(PA,PA,pB,hA,HA,HB),'mNef,ABfmJN->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(HA,hB,:,pB),t3b(PA,PA,pB,HA,HA,hB),'Mnef,ABfMJn->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(hA,hA,:,pB),t3b(PA,PA,pB,hA,HA,hB),'mnef,ABfmJn->ABeJ');

chi2A_PPvH = -D1 - D2;

X3B_ABCIJk = einsum_kg(chi2A_PPvH,t2b(:,PB,HA,hB),'ABeJ,eCIk->ABCIJk');
X3B_ABCIJk = X3B_ABCIJk - permute(X3B_ABCIJk,[1,2,3,5,4,6]);


% 
chi2A_vvvo = -0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,abfmjn->abej')...
                 -einsum_kg(sys.vB_oovv,t3b,'mnef,abfmjn->abej');

X3B_full = einsum_kg(chi2A_vvvo,t2b,'abej,ecik->abcijk');
X3B_full = X3B_full - permute(X3B_full,[1,2,3,5,4,6]);


get_error(chi2A_PPvH,chi2A_vvvo(PA,PA,:,HA))
get_error(X3B_ABCIJk,X3B_full(PA,PA,PB,HA,HA,hB))


% h2A(oooo) * t3b
D1 = 0.5*einsum_kg(H2A.oooo(HA,HA,HA,HA),t3b(PA,PA,PB,HA,HA,hB),'MNIJ,ABCMNk->ABCIJk');
D1 = D1 + einsum_kg(H2A.oooo(hA,HA,HA,HA),t3b(PA,PA,PB,hA,HA,hB),'mNIJ,ABCmNk->ABCIJk');
D1 = D1 + 0.5*einsum_kg(H2A.oooo(hA,hA,HA,HA),t3b(PA,PA,PB,hA,hA,hB),'mnIJ,ABCmnk->ABCIJk');

D1_full = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');

get_error(D1,D1_full(PA,PA,PA,HA,HA,hB))



%% ficticious t vectors
clear all
clc
close all

load h2o-631g-stretched

nfzc = 0; nfzv = 0; 
nact_h = 2; nact_p = 3; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
flag_act_scheme = 2;

sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                           nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);
                       
PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

t2a = rand(sys.Nvir_alpha,sys.Nvir_alpha,sys.Nocc_alpha,sys.Nocc_alpha);
t2b = rand(sys.Nvir_alpha,sys.Nvir_beta,sys.Nocc_alpha,sys.Nocc_beta);
t2c = rand(sys.Nvir_beta,sys.Nvir_beta,sys.Nocc_beta,sys.Nocc_beta);

t3a = rand(sys.Nvir_alpha,sys.Nvir_alpha,sys.Nvir_alpha,sys.Nocc_alpha,sys.Nocc_alpha,sys.Nocc_alpha);
t3b = rand(sys.Nvir_alpha,sys.Nvir_alpha,sys.Nvir_beta,sys.Nocc_alpha,sys.Nocc_alpha,sys.Nocc_beta);
t3c = rand(sys.Nvir_alpha,sys.Nvir_beta,sys.Nvir_beta,sys.Nocc_alpha,sys.Nocc_beta,sys.Nocc_beta);
t3d = rand(sys.Nvir_beta,sys.Nvir_beta,sys.Nvir_beta,sys.Nocc_beta,sys.Nocc_beta,sys.Nocc_beta);

t2a = asym(t2a,'a'); t2b = asym(t2b,'b'); t2c = asym(t2c,'c');
t3a = asym(t3a,'a'); t3b = asym(t3b,'b'); t3c = asym(t3c,'c'); t3d = asym(t3d,'d');

T3A.PPPHHH = t3a(PA,PA,PA,HA,HA,HA);
T3A.PPPhHH = t3a(PA,PA,PA,hA,HA,HA);
T3A.PPPhhH = t3a(PA,PA,PA,hA,hA,HA);
T3A.PPpHHH = t3a(PA,PA,pA,HA,HA,HA);
T3A.PPphHH = t3a(PA,PA,pA,hA,HA,HA);
T3A.PPphhH = t3a(PA,PA,pA,hA,hA,HA);

T3B.PPPHHH = t3b(PA,PA,PA,HA,HA,HA);
T3B.PPPhHH = t3b(PA,PA,PA,hA,HA,HA);
T3B.PPPHHh = t3b(PA,PA,PA,HA,HA,hA);
T3B.PPPhHh = t3b(PA,PA,PA,hA,HA,hA);
T3B.PPpHHH = t3b(PA,PA,pA,HA,HA,HA);
T3B.PPphHH = t3b(PA,PA,pA,hA,HA,HA);
T3B.PPpHHh = t3b(PA,PA,pA,HA,HA,hA);
T3B.PPphHh = t3b(PA,PA,pA,hA,HA,hA);


% (Vt3) diagram intermediate - chi(PPvH) (chi(ABeJ))
D1 = 0.5*einsum_kg(sys.vA_oovv(HA,HA,:,PA),T3A.PPPHHH,'MNeF,ABFMJN->ABeJ');
D1 = D1 + einsum_kg(sys.vA_oovv(hA,HA,:,PA),T3A.PPPhHH,'mNeF,ABFmJN->ABeJ');
D1 = D1 + 0.5*einsum_kg(sys.vA_oovv(hA,hA,:,PA),T3A.PPPhhH,'mneF,AFBmnJ->ABeJ');
D1 = D1 + 0.5*einsum_kg(sys.vA_oovv(HA,HA,:,pA),T3A.PPpHHH,'MNef,ABfMJN->ABeJ');
D1 = D1 + einsum_kg(sys.vA_oovv(hA,HA,:,pA),T3A.PPphHH,'mNef,ABfmJN->ABeJ');
D1 = D1 - 0.5*einsum_kg(sys.vA_oovv(hA,hA,:,pA),T3A.PPphhH,'mnef,ABfmnJ->ABeJ');

D2 = einsum_kg(sys.vB_oovv(HA,HB,:,PB),T3B.PPPHHH,'MNeF,ABFMJN->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(hA,HB,:,PB),T3B.PPPhHH,'mNeF,ABFmJN->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(HA,hB,:,PB),T3B.PPPHHh,'MneF,ABFMJn->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(hA,hB,:,PB),T3B.PPPhHh,'mneF,ABFmJn->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(HA,HB,:,pB),T3B.PPpHHH,'MNef,ABfMJN->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(hA,HB,:,pB),T3B.PPphHH,'mNef,ABfmJN->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(HA,hB,:,pB),T3B.PPpHHh,'Mnef,ABfMJn->ABeJ');
D2 = D2 + einsum_kg(sys.vB_oovv(hA,hA,:,pB),T3B.PPphHh,'mnef,ABfmJn->ABeJ');

chi2A_PPvH = -D1 - D2;

X3B_ABCIJk = einsum_kg(chi2A_PPvH,t2b(:,PB,HA,hB),'ABeJ,eCIk->ABCIJk');
X3B_ABCIJk = X3B_ABCIJk - permute(X3B_ABCIJk,[1,2,3,5,4,6]);


% 
chi2A_vvvo = -0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,abfmjn->abej')...
                 -einsum_kg(sys.vB_oovv,t3b,'mnef,abfmjn->abej');

X3B_full = einsum_kg(chi2A_vvvo,t2b,'abej,ecik->abcijk');
X3B_full = X3B_full - permute(X3B_full,[1,2,3,5,4,6]);


get_error(chi2A_PPvH,chi2A_vvvo(PA,PA,:,HA))
get_error(X3B_ABCIJk,X3B_full(PA,PA,PB,HA,HA,hB))

%%
clear all
clc
close all

nx = 5;
ny = 5;

n = 1000;
M = n/2;

lb = -500; ub = 500;
A = (ub-lb).*rand([nx,n])+lb;
B = (ub-lb).*rand([n,ny])+lb;

C = einsum_kg(A,B,'ik,kj->ij');

C2 =  einsum_kg(A(:,1:M),B(1:M,:),'iK,Kj->ij')...
     +einsum_kg(A(:,M+1:n),B(M+1:n,:),'ik,kj->ij');
 
err = abs(C-C2);
max(err(:))

