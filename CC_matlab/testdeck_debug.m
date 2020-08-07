clear all
clc
close all

format long

addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
addpath(genpath('/Users/harellab/Desktop/CC_matlab_tests/'));

%%
workpath = '/Users/karthik/Desktop/CC_matlab_tests/h2o-631g/stretched/';
[e1int, e2int, Vnuc, Norb] = load_integrals(workpath);
Nocc_a = 5;
Nocc_b = 5;
Nelec = 10;

%%
%load h2o-pvdz.mat
%load h2-pvdz-integrals.mat
%load h2o-631g-stretched
%load h2o-631g
%load rectangle-pvdz-d2h
load h2o-631g-gms

nfzc = 0; nfzv = 0; nact_h = 200; nact_p = 200;

sys_ucc = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc);
sys_cc = build_system_cc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,nfzv,nact_h,nact_p);

oa = [1:2:sys_ucc.Nelec];
ua = [sys_ucc.Nelec+1:2:2*sys_ucc.Norb]-sys_ucc.Nelec;
ob = 2:2:sys_ucc.Nelec;
ub = [sys_ucc.Nelec+2:2:2*sys_ucc.Norb]-sys_ucc.Nelec;

%% UCCSD

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-9;
ccopts.shift = 0.0;

[cc_t,Ecorr_ucc] = uccsd(sys_ucc,ccopts);

% UCCSD HBar

flag_build_3body = true;
[HBar_t] = build_ucc_HBar( cc_t, sys_ucc, flag_build_3body);
[HBar_conv, HBar_cc, t1_conv, t2_conv] = check_HBar(cc_t,sys_ucc,sys_cc);
HBar_conv{3} = HBar_cc{3};

% WE NEED TO FIGURE OUT HOW TO PERMUTE 2- AND 3-BODY HBARS WITH THE RIGHT
% SIGN AND EVERYTHING!!! UNTIL THEN, BE VERY CAREFUL AND ONLY USE THE FORM
% OF HBAR THAT YOU GET DIAGRAMMATICALLY. I.E. DO NOT USE HBARS THAT YOU GET
% FROM PERMUTING.

% Left UCCSD

lccopts.diis_size = 5;
lccopts.maxit = 500;
lccopts.tol = ccopts.tol;
lccopts.shift = 0.0;

[cc_t,lccsd_resid] = luccsd(cc_t,HBar_t,sys_ucc,lccopts);

% EOM-UCCSD

eomopts.nroot = 1;
eomopts.maxit = 100;
eomopts.tol = 1e-10;
eomopts.nvec_per_root = 1;
eomopts.max_nvec_per_root = 5;
eomopts.flag_verbose = 1;
eomopts.init_guess = 'cis';
eomopts.mult = 1;
eomopts.thresh_vec = 10*eomopts.tol;
eomopts.solver = 2;

[Rvec, omega, eom_residual, cc_t] = eomuccsd(HBar_t,cc_t,sys_ucc,eomopts);


% Left EOM-UCCSD

lccopts.maxit = 50;
lccopts.diis_size = 5;
lccopts.tol = 1e-8; 
lccopts.shift = 0;
lccopts.solver = 1;
lccopts.nroot = length(omega);

[Lvec,eom_lcc_resid,cc_t] = lefteomuccsd(omega,Rvec,HBar_t,cc_t,sys_ucc,lccopts);

%% Debug HR builds
clc

iroot = 1;

% use r ampltidues that have the correct permutational symmetries
% i.e. take t amplitudes as the values are arbitrary
r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};

% CC EOM
R1 = convert_spinint_to_spinorb({r1a,r1b},sys_ucc);
R2 = convert_spinint_to_spinorb({r2a,r2b,r2c},sys_ucc);

[HR1, HR2] = build_HBar_R(R1, R2, HBar_conv);

% UCC EOM
[X1A] = build_HR_1A(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys_ucc);
[X1B] = build_HR_1B(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys_ucc);
[X2A] = build_HR_2A(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys_ucc);
[X2B] = build_HR_2B(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys_ucc);
[X2C] = build_HR_2C(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys_ucc);

X1 = convert_spinint_to_spinorb({X1A,X1B},sys_ucc);
X2 = convert_spinint_to_spinorb({X2A,X2B,X2C},sys_ucc);

% can also debug specific diagrams
% get_error(Y(ua,ua,oa,oa),Y2A)
% get_error(Y(ua,ub,oa,ob),Y2B)
% get_error(Y(ub,ub,ob,ob),Y2C)

%
fprintf('Error in UCC HR on singles = %4.10f\n',get_error(X1,HR1))
fprintf('Error in UCC HR on doubles = %4.10f\n',get_error(X2,HR2))

%% Debug LH builds
clc

iroot = 1;

if iroot == 0
    flag_ground = true; 
    iroot2 = 1;
    omega = 0.0;
else
    flag_ground = false;
    iroot2 = iroot+1;
    omega(iroot2-1);
end

shift = 0.0;

l1a = cc_t.l1a{iroot2}; l1b = cc_t.l1b{iroot2};
l2a = cc_t.l2a{iroot2}; l2b = cc_t.l2b{iroot2}; l2c = cc_t.l2c{iroot2};

% CC left
[L1] = convert_spinint_to_spinorb({l1a,l1b},sys_ucc);
[L2] = convert_spinint_to_spinorb({l2a,l2b,l2c},sys_ucc);
%
flag_jacobi = 1;
[LH1, LH2] = build_L_HBar(L1, L2, HBar_conv, flag_ground, flag_jacobi);

% UCC left
flag_jacobi = true;
[X1A] = build_LH_1A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys_ucc,flag_jacobi,flag_ground);
[X1B] = build_LH_1B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys_ucc,flag_jacobi,flag_ground);
[X2A] = build_LH_2A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys_ucc,flag_jacobi,flag_ground);
[X2B] = build_LH_2B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys_ucc,flag_jacobi,flag_ground);
[X2C] = build_LH_2C(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys_ucc,flag_jacobi,flag_ground);

X1 = convert_spinint_to_spinorb({X1A,X1B},sys_ucc);
X2 = convert_spinint_to_spinorb({X2A,X2B,X2C},sys_ucc);

%
fprintf('Error in UCC LH on singles = %4.10f\n',get_error(X1,LH1))
fprintf('Error in UCC LH on doubles = %4.10f\n',get_error(X2,LH2))

%% Debug MM23 build
clc

t1 = convert_spinint_to_spinorb({cc_t.t1a,cc_t.t1b},sys_ucc);
t2 = convert_spinint_to_spinorb({cc_t.t2a,cc_t.t2b,cc_t.t2c},sys_ucc);

[MM23] = build_MM23(t1,t2,HBar_conv);

MM23A_conv = MM23(ua,ua,ua,oa,oa,oa);
MM23B_conv = MM23(ua,ua,ub,oa,oa,ob);
MM23C_conv = MM23(ua,ub,ub,oa,ob,ob);
MM23D_conv = MM23(ub,ub,ub,ob,ob,ob);

[MM23A] = build_MM23A(cc_t,HBar_t,sys_ucc);
[MM23B] = build_MM23B(cc_t,HBar_t,sys_ucc);
[MM23C] = build_MM23C(cc_t,HBar_t,sys_ucc);
[MM23D] = build_MM23D(cc_t,HBar_t,sys_ucc);

fprintf('Error in MM23A = %4.10f\n',get_error(MM23A,MM23A_conv))
fprintf('Error in MM23B = %4.10f\n',get_error(MM23B,MM23B_conv))
fprintf('Error in MM23C = %4.10f\n',get_error(MM23C,MM23C_conv))
fprintf('Error in MM23D = %4.10f\n',get_error(MM23D,MM23D_conv))

%% Debug approximate L3 build
clc

iroot = 1;

if iroot == 0
    flag_ground = true; 
    iroot2 = 1;
    omega = 0.0;
else
    flag_ground = false;
    iroot2 = iroot+1;
    omega(iroot2-1);
end

% CC left
[L1] = convert_spinint_to_spinorb({cc_t.l1a{iroot2},cc_t.l1b{iroot2}},sys_ucc);
[L2] = convert_spinint_to_spinorb({cc_t.l2a{iroot2},cc_t.l2b{iroot2},cc_t.l2c{iroot2}},sys_ucc);

[L3] = build_L3_approx(L1,L2,HBar_conv);
L3A_conv = L3(oa,oa,oa,ua,ua,ua);
L3B_conv = L3(oa,oa,ob,ua,ua,ub);
L3C_conv = L3(oa,ob,ob,ua,ub,ub);
L3D_conv = L3(ob,ob,ob,ub,ub,ub);

[L3A] = build_L3A_approx(cc_t,HBar_t,sys_ucc,iroot);
[L3B] = build_L3B_approx(cc_t,HBar_t,sys_ucc,iroot);
[L3C] = build_L3C_approx(cc_t,HBar_t,sys_ucc,iroot);
[L3D] = build_L3D_approx(cc_t,HBar_t,sys_ucc,iroot);

fprintf('Error in L3A = %4.10f\n',get_error(L3A,L3A_conv))
fprintf('Error in L3B = %4.10f\n',get_error(L3B,L3B_conv))
fprintf('Error in L3C = %4.10f\n',get_error(L3C,L3C_conv))
fprintf('Error in L3D = %4.10f\n',get_error(L3D,L3D_conv))


%% Debug new 3-body HBars for EOM-MM23
clc
H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3C = HBar_t.H3C; H3D = HBar_t.H3D;
H3 = HBar_cc{3};

fprintf('=============H3(vvvoov)=============\n')
fprintf('Error in H3A(vvvoov) = %4.10f\n',get_error(H3A.vvvoov,H3{2,2,2,1,1,2}(ua,ua,ua,oa,oa,ua)))
fprintf('Error in H3B(vvvoov) = %4.10f\n',get_error(H3B.vvvoov,H3{2,2,2,1,1,2}(ua,ua,ub,oa,oa,ub)))
fprintf('Error in H3B(vvvovo) = %4.10f\n',get_error(H3B.vvvovo,H3{2,2,2,1,2,1}(ua,ua,ub,oa,ua,ob)))
fprintf('Error in H3C(vvvoov) = %4.10f\n',get_error(H3C.vvvoov,H3{2,2,2,1,1,2}(ua,ub,ub,oa,ob,ub)))
fprintf('Error in H3C(vvvvoo) = %4.10f\n',get_error(H3C.vvvvoo,H3{2,2,2,2,1,1}(ua,ub,ub,ua,ob,ob)))
fprintf('Error in H3D(vvvoov) = %4.10f\n',get_error(H3D.vvvoov,H3{2,2,2,1,1,2}(ub,ub,ub,ob,ob,ub)))
fprintf('====================================\n')
fprintf('\n')
fprintf('=============H3(vvoooo)=============\n')
fprintf('Error in H3A(vvoooo) = %4.10f\n',get_error(H3A.vvoooo,H3{2,2,1,1,1,1}(ua,ua,oa,oa,oa,oa)))
fprintf('Error in H3B(vvoooo) = %4.10f\n',get_error(H3B.vvoooo,H3{2,2,1,1,1,1}(ua,ua,ob,oa,oa,ob)))
fprintf('Error in H3B(vovooo) = %4.10f\n',get_error(H3B.vovooo,H3{2,1,2,1,1,1}(ua,oa,ub,oa,oa,ob)))
fprintf('Error in H3C(ovvooo) = %4.10f\n',get_error(H3C.ovvooo,H3{1,2,2,1,1,1}(oa,ub,ub,oa,ob,ob)))
fprintf('Error in H3C(vvoooo) = %4.10f\n',get_error(H3C.vvoooo,H3{2,2,1,1,1,1}(ua,ub,ob,oa,ob,ob)))
fprintf('Error in H3D(vvoooo) = %4.10f\n',get_error(H3D.vvoooo,H3{2,2,1,1,1,1}(ub,ub,ob,ob,ob,ob)))
fprintf('====================================\n')
fprintf('\n')
fprintf('=============H3(vvvvvo)=============\n')
fprintf('Error in H3A(vvvvvo) = %4.10f\n',get_error(H3A.vvvvvo,H3{2,2,2,2,2,1}(ua,ua,ua,ua,ua,oa)))
fprintf('Error in H3B(vvvovv) = %4.10f\n',get_error(H3B.vvvovv,H3{2,2,2,1,2,2}(ua,ua,ub,oa,ua,ub)))
fprintf('Error in H3B(vvvvvo) = %4.10f\n',get_error(H3B.vvvvvo,H3{2,2,2,2,2,1}(ua,ua,ub,ua,ua,ob)))
fprintf('Error in H3C(vvvvvo) = %4.10f\n',get_error(H3C.vvvvvo,H3{2,2,2,2,2,1}(ua,ub,ub,ua,ub,ob)))
fprintf('Error in H3C(vvvovv) = %4.10f\n',get_error(H3C.vvvovv,H3{2,2,2,1,2,2}(ua,ub,ub,oa,ub,ub)))
fprintf('Error in H3D(vvvvvo) = %4.10f\n',get_error(H3D.vvvvvo,H3{2,2,2,2,2,1}(ub,ub,ub,ub,ub,ob)))
fprintf('====================================\n')
fprintf('\n')
fprintf('=============H3(oovooo)=============\n')
fprintf('Error in H3A(oovooo) = %4.10f\n',get_error(H3A.oovooo,H3{1,1,2,1,1,1}(oa,oa,ua,oa,oa,oa)))
fprintf('Error in H3B(oovooo) = %4.10f\n',get_error(H3B.oovooo,H3{1,1,2,1,1,1}(oa,oa,ub,oa,oa,ob)))
fprintf('Error in H3B(vooooo) = %4.10f\n',get_error(H3B.vooooo,H3{2,1,1,1,1,1}(ua,oa,ob,oa,oa,ob)))
fprintf('Error in H3C(oovooo) = %4.10f\n',get_error(H3C.oovooo,H3{1,1,2,1,1,1}(oa,ob,ub,oa,ob,ob)))
fprintf('Error in H3C(vooooo) = %4.10f\n',get_error(H3C.vooooo,H3{2,1,1,1,1,1}(ua,ob,ob,oa,ob,ob)))
fprintf('Error in H3D(oovooo) = %4.10f\n',get_error(H3D.oovooo,H3{1,1,2,1,1,1}(ob,ob,ub,ob,ob,ob)))
fprintf('====================================\n')

%% Debug EOM-MM23 builds
clc

iroot = 1;

r0 = cc_t.r0(iroot);

r1a = cc_t.r1a{iroot}; r1b = cc_t.r1b{iroot};
r2a = cc_t.r2a{iroot}; r2b = cc_t.r2b{iroot}; r2c = cc_t.r2c{iroot};

r1 = convert_spinint_to_spinorb({r1a,r1b},sys_ucc);
r2 = convert_spinint_to_spinorb({r2a,r2b,r2c},sys_ucc);

t1 = convert_spinint_to_spinorb({cc_t.t1a,cc_t.t1b},sys_ucc);
t2= convert_spinint_to_spinorb({cc_t.t2a,cc_t.t2b,cc_t.t2c},sys_ucc);

[EOMMM23] = build_EOMMM23(t1,t2,r1,r2,HBar_conv);

EOMMM23A_conv = EOMMM23(ua,ua,ua,oa,oa,oa);
EOMMM23B_conv = EOMMM23(ua,ua,ub,oa,oa,ob);
EOMMM23C_conv = EOMMM23(ua,ub,ub,oa,ob,ob);
EOMMM23D_conv = EOMMM23(ub,ub,ub,ob,ob,ob);

[EOMMM23A] = build_EOMMM23A(cc_t,HBar_t,iroot);
[EOMMM23B] = build_EOMMM23B(cc_t,HBar_t,iroot);
[EOMMM23C] = build_EOMMM23C(cc_t,HBar_t,iroot);
[EOMMM23D] = build_EOMMM23D(cc_t,HBar_t,iroot);


fprintf('Error in EOMMM23A = %4.10f\n',get_error(EOMMM23A,EOMMM23A_conv))
fprintf('Error in EOMMM23B = %4.10f\n',get_error(EOMMM23B,EOMMM23B_conv))
fprintf('Error in EOMMM23C = %4.10f\n',get_error(EOMMM23C,EOMMM23C_conv))
fprintf('Error in EOMMM23D = %4.10f\n',get_error(EOMMM23D,EOMMM23D_conv))


%%
sys = sys_cc;

D1 =        einsum_kg(sys.Vvvvv,t2,'acfe,fbij->abcije')...
           -einsum_kg(sys.Vvvvv,t2,'bcfe,faij->abcije')...
           -einsum_kg(sys.Vvvvv,t2,'abfe,fcij->abcije');
       
D2 = einsum_kg(einsum_kg(sys.Voovo,t1,'mnei,cm->cnei'),t2,'cnei,abnj->abcije')...
     -einsum_kg(einsum_kg(sys.Voovo,t1,'mnei,am->anei'),t2,'anei,cbnj->abcije')... % (ac)
     -einsum_kg(einsum_kg(sys.Voovo,t1,'mnei,bm->bnei'),t2,'bnei,acnj->abcije')... % (bc)
     -einsum_kg(einsum_kg(sys.Voovo,t1,'mnej,cm->cnej'),t2,'cnej,abni->abcije')... % (ij)
     +einsum_kg(einsum_kg(sys.Voovo,t1,'mnei,am->anej'),t2,'anej,cbni->abcije')... % (ij)(ac)
     +einsum_kg(einsum_kg(sys.Voovo,t1,'mnej,bm->bnej'),t2,'bnej,acni->abcije');   % (ij)(bc)

D3 = -einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fi->amie'),t2,'amie,bcjm->abcije')...
     +einsum_kg(einsum_kg(sys.Vvovv,t1,'bmfe,fi->bmie'),t2,'bmie,acjm->abcije')... % (ab)
     +einsum_kg(einsum_kg(sys.Vvovv,t1,'cmfe,fi->cmie'),t2,'cmie,bajm->abcije')... % (ac)
     +einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,fj->amje'),t2,'amje,bcim->abcije')... % (ij)
     -einsum_kg(einsum_kg(sys.Vvovv,t1,'bmfe,fj->bmje'),t2,'bmje,acim->abcije')... % (ij)(ab)
     -einsum_kg(einsum_kg(sys.Vvovv,t1,'cmfe,fj->cmje'),t2,'cmje,baim->abcije');   % (ij)(ac)
    
D4 = -einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,cm->acfe'),t2,'acfe,fbij->abcije')...
     +einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,cm->acfe'),t2,'acfe,fbij->abcije')...
     +einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,cm->acfe'),t2,'acfe,fbij->abcije')...
     +einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,cm->acfe'),t2,'acfe,fbij->abcije')...
     -einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,cm->acfe'),t2,'acfe,fbij->abcije')...
     -einsum_kg(einsum_kg(sys.Vvovv,t1,'amfe,cm->acfe'),t2,'acfe,fbij->abcije');



