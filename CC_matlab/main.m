clear all
clc
close all

format long

%addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab/'));
%addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/Davidson_matlab_test'));
%addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/GDav'));
%addpath(genpath('/Users/karthik/Dropbox/Davidson_methods/davidson_cheb_test/Bchdav'));

% addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2o-pvdz/'));
% addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2o-631g/'));
% addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/h2-pvdz/'));
% addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/'));

addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests'));

%addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));
%addpath(genpath('/Users/harellab/Desktop/CC_matlab_tests/'));

%%
workpath = '/Users/karthik/Desktop/CC_matlab_tests/square-d2h-pvdz/';
[e1int, e2int, Vnuc, Norb] = load_integrals(workpath);
Nocc_a = 14;
Nocc_b = 14;
Nelec = 28;

%%
%load h2o-pvdz.mat
%load h2-pvdz-integrals.mat
%load h2o-631g-stretched
%load h2o-631g
%load rectangle-pvdz-d2h % NEEDS NFZC = 4    
%load square-pvdz-d2h % NEEDS NFZC = 4
%load h2o-pvdz-gms.mat
%load h2o-631g-gms.mat
load f2-2.0-pvdz.mat % NEEDS NFZC = 2
%load f2-1.0-pvdz.mat % NEEDS NFZC = 2

%%
nfzc = 2; nfzv = 0; nact = 100;
nact_h_alpha = nact; nact_h_beta = nact; nact_p_alpha = nact; nact_p_beta = nact;
flag_act_scheme = 0;
%nact_h = 200; nact_p = 200;
sys_ucc = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                         nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);
%sys_cc = build_system_cc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,nfzv,nact_h,nact_p);
%clear ans e1int e2int


%% CIS

cisopts.nroot = 10;
cisopts.maxit = 50;
cisopts.tol = 1e-8;
cisopts.nvec_per_root = 1;
cisopts.max_nvec_per_root = 5;
cisopts.flag_verbose = 1;
cisopts.thresh_vec = 1e-3;

[C, omega_cis, c0] = cis(sys_ucc,cisopts);

%% UCCSD

ccopts.diis_size = 5;
ccopts.maxit = 200;
ccopts.tol = 1e-8;
ccopts.shift = 0;
[cc_t,Ecorr_uccsd] = uccsd(sys_ucc,ccopts);

%% UCCSDT


% H2O / 631G (C1) - stretched 
% E(Ref)=  -75.7199662951565
% E(Cor)= -0.221361597785899
% E(CCSDT)=  -75.9413278929424

% F2 / cc-pVDZ (D2h) - R/Re = 2.0 
% E(Ref)= -198.420096281519
% E(Cor)= -0.638105013177
% E(CCSDT)= -199.058201294696

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-10;
ccopts.shift = 0;

[cc_t,Ecorr_uccsdt] = uccsdt(sys_ucc,ccopts);

%% Active space CCSDt-3 variant

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-10;
ccopts.shift = 0;

[cc_t,Ecorr_uccsdt] = uccsdt3(sys_ucc,ccopts);

%% UCCSD HBar

flag_3body = false;
[HBar_t] = build_ucc_HBar( cc_t, sys_ucc, flag_3body);

%% Left UCCSD

lccopts.diis_size = 5;
lccopts.maxit = 500;
lccopts.tol = ccopts.tol;
lccopts.shift = 0.0;

[cc_t,lccsd_resid] = luccsd(cc_t,HBar_t,sys_ucc,lccopts);


%% EOM-UCCSD

eomopts.nroot = 5;
eomopts.maxit = 100;
eomopts.tol = 1e-6;
eomopts.nvec_per_root = 1;
eomopts.max_nvec_per_root = 5;
eomopts.flag_verbose = 1;
eomopts.init_guess = 'cis';
eomopts.mult = 1;
eomopts.thresh_vec = 1e-3;
eomopts.solver = 2;

[Rvec, omega, eom_residual, cc_t] = eomuccsd(HBar_t,cc_t,sys_ucc,eomopts);


%% Left EOM-UCCSD

lccopts.maxit = 50;
lccopts.diis_size = 5;
lccopts.tol = 1e-8; 
lccopts.shift = 0;
lccopts.solver = 1;
lccopts.nroot = length(omega);

[Lvec,eom_lcc_resid,cc_t] = lefteomuccsd(omega,Rvec,HBar_t,cc_t,sys_ucc,lccopts);

%% CR-CC(2,3) 

[Ecrcc23A,Ecrcc23B,Ecrcc23C,Ecrcc23D] = crcc23_wrap(cc_t,HBar_t,sys_ucc);

%%
%%%%%%%%%%%%%%%%%%%% Spinorbital CC Codes %%%%%%%%%%%%%%%%%%%%
%
%% CCSD

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-10;
ccopts.shift = 0.0;

[t1,t2,Ecorr_ccsd] = ccsd(sys_cc,ccopts);

%% CCSDT

ccopts.diis_size = 5;
ccopts.maxit = 100;
ccopts.tol = 1e-10;
ccopts.shift = 0;

% CCSDT errors on order of 10^-5... be VERY careful checking all antisymmetrizers with einsum!
[t1,t2,t3,Ecorr_ccsdt] = ccsdt(sys_cc,ccopts);
%[t1,t2,t3,Ecorr_ccsdt3] = ccsdt3(sys_cc,ccopts);


%% Build HBar

%[HBar] = build_HBar(t1,t2,sys_cc); 
[HBar] = build_HBar_debug(t1,t2,sys_cc);

%% Left-CCSD 

lccopts.diis_size = 5;
lccopts.maxit = 200;
lccopts.tol = 1e-9;
lccopts.shift = 0.0;

[lambda1,lambda2,lcc_resid] = lccsd(t1,t2,HBar,sys_cc,lccopts);

%% CR-CC(2,3)

[crcc23A, crcc23B, crcc23C, crcc23D] = crcc23(t1,t2,lambda1,lambda2,HBar,sys_cc);
%[Ecorr_crcc23A,Ecorr_crcc23B,Ecorr_crcc23C,Ecorr_crcc23D] = crcc23_opt(t1,t2,lambda1,lambda2,HBar,sys_cc);

%% EOM-CCSD

eomopts.nroot = 5;
eomopts.maxit = 100;
eomopts.tol = 1e-6;
eomopts.nvec_per_root = 1;
eomopts.max_nvec_per_root = 5;
eomopts.flag_verbose = 1;
eomopts.init_guess = 'cis';
eomopts.thresh_vec = 10*eomopts.tol;


% it's possible there's an error in EOMCCSD sigma function. R vector for
% first singlet state of H2 are close to Jun's but not exactly identical...
[Rvec, omega, r0, eom_residual] = eomccsd(HBar,t1,t2,sys_cc,sys_ucc,eomopts);

% Using iterative diagonalization solvers of Yunkai Zhou
% these also encounter subspace collapse and 0 eigenvalue problems!

% eomopts.tol = 1e-4;
% eomopts.Matsymm = 1;
% eomopts.vmax = max(2*nroot, 20);
% eomopts.itmax = 100;
% eomopts.nkeep = nroot;
% eomopts.v0 = rand(Nocc*Nunocc+Nocc^2*Nunocc^2,1);
% eomopts.displ = 2;
% eomopts.les_solver = 31; % 71 jacobi-davidson minres does not work!

% HRvec = @(x) HR_matvec(x,HBar);
% [eval, V, nconv, history] = gdav(HRvec, Nocc*Nunocc+Nocc^2*Nunocc^2, eomopts.nroot, 'SM',eomopts);
% [eval, V, nconv, history] = bchdav('HR_matvec', Nocc*Nunocc+Nocc^2*Nunocc^2, eomopts.nroot);

%% Left-EOMCCSD 

lccopts.maxit = 80;
lccopts.diis_size = 5;
lccopts.tol = eomopts.tol; % can only converge Left-EOM to the same tolerance as we converged right EOM
lccopts.shift = 0.00;

[Lvec,omega_lcc,eom_lcc_resid] = lefteomccsd(omega,Rvec,HBar,t1,t2,sys_cc,lccopts);

%% CR-EOMCC(2,3)

% don't run with spinorbital solvers. Use spin-integrated CCSD through
% EOM and left and convert to spinorbital

iroot = 2;

l1a = reshape(Lvec(sys_ucc.posv{1},iroot),sys_ucc.size{1});
l1b = reshape(Lvec(sys_ucc.posv{2},iroot),sys_ucc.size{2});
l2a = reshape(Lvec(sys_ucc.posv{3},iroot),sys_ucc.size{3});
l2b = reshape(Lvec(sys_ucc.posv{4},iroot),sys_ucc.size{4});
l2c = reshape(Lvec(sys_ucc.posv{5},iroot),sys_ucc.size{5});

r1a = reshape(Rvec(sys_ucc.posv{1},iroot),sys_ucc.size{1});
r1b = reshape(Rvec(sys_ucc.posv{2},iroot),sys_ucc.size{2});
r2a = reshape(Rvec(sys_ucc.posv{3},iroot),sys_ucc.size{3});
r2b = reshape(Rvec(sys_ucc.posv{4},iroot),sys_ucc.size{4});
r2c = reshape(Rvec(sys_ucc.posv{5},iroot),sys_ucc.size{5});

t1 = convert_spinint_to_spinorb({cc_t.t1a,cc_t.t1b},sys_ucc);
t2 = convert_spinint_to_spinorb({cc_t.t2a,cc_t.t2b,cc_t.t2c},sys_ucc);
r1 = convert_spinint_to_spinorb({r1a,r1b},sys_ucc);
r2 = convert_spinint_to_spinorb({r2a,r2b,r2c},sys_ucc);
l1 = convert_spinint_to_spinorb({l1a,l1b},sys_ucc);
l2 = convert_spinint_to_spinorb({l2a,l2b,l2c},sys_ucc);

flag_build_3body = true;
[HBar_t] = build_ucc_HBar( cc_t, sys_ucc, flag_build_3body);
[HBar_conv, HBar_cc, ~, ~] = check_HBar(cc_t,sys_ucc,sys_cc);
HBar_conv{3} = HBar_cc{3};

[omega_crcc23A,omega_crcc23B,omega_crcc23C,omega_crcc23D,LM3] = creomcc23(t1,t2,r1,r2,l1,l2,r0(iroot),omega(iroot),HBar_conv,sys_cc);



%% Checking Jun's L vectors vs. spinint
clc

fid = fopen('L-CCSDt');
L = fread(fid,Inf,'double');
fclose(fid)

iroot = [0, 2, 5, 6];
%iroot_jun = [1, 2, 3, 4];
% iroot_jun = [1, 2];
% iroot = [0, 2];


iroot_jun = [1];
for idx = 1:length(iroot_jun)
    
    J = iroot_jun(idx);
    j = iroot(idx);

    L_jun = L((J-1)*sys_ucc.triples_dim+1:J*sys_ucc.triples_dim);
    l1a_jun = reshape(L_jun(sys_ucc.posv{1}),sys_ucc.size{1});
    l1b_jun = reshape(L_jun(sys_ucc.posv{2}),sys_ucc.size{2});
    l2a_jun = reshape(L_jun(sys_ucc.posv{3}),sys_ucc.size{3});
    l2b_jun = reshape(L_jun(sys_ucc.posv{4}),sys_ucc.size{4});
    l2c_jun = reshape(L_jun(sys_ucc.posv{5}),sys_ucc.size{5});

    l1_jun = convert_spinint_to_spinorb({l1a_jun,l1b_jun},sys_ucc);
    l2_jun = convert_spinint_to_spinorb({l2a_jun,l2b_jun,l2c_jun},sys_ucc);

    if j == 0
        l1a = cc_t.l1a{1};
        l1b = cc_t.l1b{1};
        l2a = cc_t.l2a{1};
        l2b = cc_t.l2b{1};
        l2c = cc_t.l2c{1};
    else
        l1a = cc_t.l1a{j+1};
        l1b = cc_t.l1b{j+1};
        l2a = cc_t.l2a{j+1};
        l2b = cc_t.l2b{j+1};
        l2c = cc_t.l2c{j+1};
    end

    l1 = convert_spinint_to_spinorb({l1a,l1b},sys_ucc);
    l2 = convert_spinint_to_spinorb({l2a,l2b,l2c},sys_ucc);
    
    if j ~= 0
        fprintf('\nError on root %d:  VEE = %4.8f\n',idx,omega(j))    

        err = get_error(abs(l1a),abs(l1a_jun));
        percerr = norm(abs(l1a(:))-abs(l1a_jun(:)))/norm(abs(l1a_jun(:)))*100;
        max_err = max(abs(abs(l1a(:))-abs(l1a_jun(:))));
        fprintf('Error(L1A) = %4.12f (%4.8f%%)\n',err,percerr)
        fprintf('      Max Error = %4.12f\n',max_err)

        err = get_error(abs(l1b),abs(l1b_jun));
        percerr = norm(abs(l1b(:))-abs(l1b_jun(:)))/norm(abs(l1b_jun(:)))*100;
        max_err = max(abs(abs(l1b(:))-abs(l1b_jun(:))));
        fprintf('Error(L1B) = %4.12f (%4.8f%%)\n',err,percerr)
        fprintf('      Max Error = %4.12f\n',max_err)

        err = get_error(abs(l2a),abs(l2a_jun));
        percerr = norm(abs(l2a(:))-abs(l2a_jun(:)))/norm(abs(l2a_jun(:)))*100;
        max_err = max(abs(l2a(:))-abs(l2a_jun(:)));
        fprintf('Error(L2A) = %4.12f (%4.8f%%)\n',err,percerr)
        fprintf('      Max Error = %4.12f\n',max_err)

        err = get_error(abs(l2b),abs(l2b_jun));
        percerr = norm(abs(l2b(:))-abs(l2b_jun(:)))/norm(abs(l2b_jun(:)))*100;
        max_err = max(abs(l2b(:))-abs(l2b_jun(:)));
        fprintf('Error(L2B) = %4.12f (%4.8f%%)\n',err,percerr)
        fprintf('      Max Error = %4.12f\n',max_err)

        err = get_error(abs(l2c),abs(l2c_jun));
        percerr = norm(abs(l2c(:))-abs(l2c_jun(:)))/norm(abs(l2c_jun(:)))*100;
        max_err = max(abs(l2c(:))-abs(l2c_jun(:)));
        fprintf('Error(L2C) = %4.12f (%4.8f%%)\n',err,percerr)
        fprintf('      Max Error = %4.12f\n',max_err)

    else 
        fprintf('Error on ground state: E0 = %4.12f\n',sys_ucc.Escf+Ecorr_uccsd)
       %fprintf('Error(L1) = %4.12f\n',get_error(l1,l1_jun))
       %fprintf('Error(L2) = %4.12f\n',get_error(l2,l2_jun))
        fprintf('Error(L1A) = %4.12f\n',get_error(l1a,l1a_jun))
        fprintf('Error(L1B) = %4.12f\n',get_error(l1b,l1b_jun))
        fprintf('Error(L2A) = %4.12f\n',get_error(l2a,l2a_jun))
        fprintf('Error(L2B) = %4.12f\n',get_error(l2b,l2b_jun))
        fprintf('Error(L2C) = %4.12f\n',get_error(l2c,l2c_jun))
    end

    

end




%%
clc

iroot = 1;

%pr0 = @(i,j,k,a,b,c) print_moments(MM23A,L3A,sys_ucc,cc_t,HBar_t,0.0,i,j,k,a,b,c);
prA = @(i,j,k,a,b,c) print_moments(EOMMM23A,EOML3A,sys_ucc,cc_t,HBar_t,omega(iroot),i,j,k,a,b,c,'A');
prB = @(i,j,k,a,b,c) print_moments(EOMMM23B,EOML3B,sys_ucc,cc_t,HBar_t,omega(iroot),i,j,k,a,b,c,'B');
prC = @(i,j,k,a,b,c) print_moments(M3C_v2,L3C_v2,sys_ucc,cc_t,HBar_t,omega(iroot),i,j,k,a,b,c,'C');
prD = @(i,j,k,a,b,c) print_moments(EOMMM23D,EOML3D,sys_ucc,cc_t,HBar_t,omega(iroot),i,j,k,a,b,c,'D');


% excited_state A
% prA(2,3,4,6,7,11)

% excited state B
% prB(1,2,1,6,9,10)
% prB(1,2,1,6,10,6)
% prB(1,2,1,6,10,9)
% prB(1,2,1,6,10,11)
% prB(1,2,1,6,10,13)

% excited state C
% prC(5,4,5,10,9,11)
% prC(5,4,5,10,7,8)
% prC(5,4,5,13,10,11)

% excited state D
% prD(3,4,5,11,12,13)
% prD(3,4,5,9,12,13)
% prD(1,2,3,8,10,11)
% prD(1,2,3,10,11,12)


% [Lsort,idx] = sort(abs(EOML3A(:)),'descend');
% N = 10;
% ct = 1;
% for p = 1:36:36*N;
%     [i,j,k,a,b,c] = ind2sub(size(EOML3A),idx(p));
%     A = a + sys_ucc.Nocc_alpha;
%     B = b + sys_ucc.Nocc_alpha;
%     C = c + sys_ucc.Nocc_alpha;
%     pr2(i,j,k,A,B,C)
%     ct = ct + 1;
% end

% ground state
% pr1(1,2,3,6,7,9)
% pr1(1,2,3,6,7,11)
% pr1(1,2,3,6,7,13)
% pr1(1,2,3,6,8,9)
% pr1(1,2,3,6,8,11)
% pr1(1,2,3,6,8,13)
% pr1(1,2,3,6,9,12)
%%

fidMA = fopen('~/Desktop/CC_matlab_tests/h2o-631g/stretched/M3A');
fidMB = fopen('~/Desktop/CC_matlab_tests/h2o-631g/stretched/M3B');
fidMC = fopen('~/Desktop/CC_matlab_tests/h2o-631g/stretched/M3C');
fidMD = fopen('~/Desktop/CC_matlab_tests/h2o-631g/stretched/M3D');

[M3A] = build_MM23A(cc_t,HBar_t,sys_ucc);
[M3B] = build_MM23B(cc_t,HBar_t,sys_ucc);
[M3C] = build_MM23C(cc_t,HBar_t,sys_ucc);
[M3D] = build_MM23D(cc_t,HBar_t,sys_ucc);

M3A_lin = zeros(get_number_unique('3A',sys_ucc),1);
M3A_jun = zeros(get_number_unique('3A',sys_ucc),1);
M3B_lin = zeros(get_number_unique('3B',sys_ucc),1);
M3B_jun = zeros(get_number_unique('3B',sys_ucc),1);
M3C_lin = zeros(get_number_unique('3C',sys_ucc),1);
M3C_jun = zeros(get_number_unique('3C',sys_ucc),1);
M3D_lin = zeros(get_number_unique('3D',sys_ucc),1);
M3D_jun = zeros(get_number_unique('3D',sys_ucc),1);

% we're not accessing M3B and M3C correctly so there's a mismatch between
% mine and Jun's..

ct = 1;
for i = 1:sys_ucc.Nocc_alpha
    for j = i+1:sys_ucc.Nocc_alpha
        for k = j+1:sys_ucc.Nocc_alpha
            for a = 1:sys_ucc.Nvir_alpha
                for b = a+1:sys_ucc.Nvir_alpha
                    for c = b+1:sys_ucc.Nvir_alpha

                        tlineMA = fgetl(fidMA);

                        M3A_jun(ct) = str2double(tlineMA);
                        M3A_lin(ct) = M3A(a,b,c,i,j,k);
                        
%                         if abs(M3A_jun(ct)) > 1e-6
%                             fprintf('M3A(%d) = %4.8f\n',ct,M3A_lin(ct))
%                             fprintf('M3A_jun(%d) = %4.8f\n',ct,M3A_jun(ct))   
%                         end
                        
                        ct = ct + 1;
                    end
                end
            end
         end
     end
end

fprintf('Error in M3A: %4.12f\n',get_error(M3A_lin,M3A_jun))

ct = 1;
for i = 1:sys_ucc.Nocc_alpha
    for j = i+1:sys_ucc.Nocc_alpha
        for k = 1:sys_ucc.Nocc_beta
            for a = 1:sys_ucc.Nvir_alpha
                for b = a+1:sys_ucc.Nvir_alpha
                    for c = 1:sys_ucc.Nvir_beta

                        tlineMB = fgetl(fidMB); 

                        M3B_jun(ct) = str2double(tlineMB);
                        M3B_lin(ct) = M3B(a,b,c,i,j,k);
                        
%                         if abs(M3B_jun(ct)) > 1e-6
%                             fprintf('M3B(%d) = %4.8f\n',ct,M3B_lin(ct))
%                             fprintf('M3B_jun(%d) = %4.8f\n',ct,M3B_jun(ct))   
%                         end
                        
                        ct = ct + 1;
                    end
                end
            end
         end
     end
end

fprintf('Error in M3B: %4.12f\n',get_error(M3B_lin,M3B_jun))

ct = 1;
for i = 1:sys_ucc.Nocc_alpha
    for j = 1:sys_ucc.Nocc_beta
        for k = i+1:sys_ucc.Nocc_beta
            for a = 1:sys_ucc.Nvir_alpha
                for b = 1:sys_ucc.Nvir_beta
                    for c = b+1:sys_ucc.Nvir_beta

                        tlineMC = fgetl(fidMC); 

                        M3C_jun(ct) = str2double(tlineMC);
                        M3C_lin(ct) = M3C(a,b,c,i,j,k);
                        
                        if abs(M3C_jun(ct)) > 1e-6
                            fprintf('M3C(%d) = %4.8f\n',ct,M3C_lin(ct))
                            fprintf('M3C_jun(%d) = %4.8f\n',ct,M3C_jun(ct))   
                        end
                        
                        ct = ct + 1;
                    end
                end
            end
         end
     end
end

fprintf('Error in M3C: %4.12f\n',get_error(M3C_lin,M3C_jun))

ct = 1;
for i = 1:sys_ucc.Nocc_beta
    for j = i+1:sys_ucc.Nocc_beta
        for k = j+1:sys_ucc.Nocc_beta
            for a = 1:sys_ucc.Nvir_beta
                for b = a+1:sys_ucc.Nvir_beta
                    for c = b+1:sys_ucc.Nvir_beta

                        tlineMD = fgetl(fidMD); 

                        M3D_jun(ct) = str2double(tlineMD);
                        M3D_lin(ct) = M3D(a,b,c,i,j,k);
                        
%                         if abs(M3D_jun(ct)) > 1e-6
%                             fprintf('M3D(%d) = %4.8f\n',ct,M3D_lin(ct))
%                             fprintf('M3D_jun(%d) = %4.8f\n',ct,M3D_jun(ct))   
%                         end
                        
                        ct = ct + 1;
                    end
                end
            end
         end
     end
end

fprintf('Error in M3D: %4.12f\n',get_error(M3D_lin,M3D_jun))
%%
ct = 1;
fidD = fopen('~/Desktop/CC_matlab_tests/h2o-631g/stretched/L3D'); fidA = fopen('L3A'); fidB = fopen('L3B'); fidC = fopen('L3C');
fidDD = fopen('D3D'); fidDA = fopen('D3A'); fidDB = fopen('D3B'); fidDC = fopen('D3C');
fidMD = fopen('M3D'); fidMA = fopen('M3A'); fidMB = fopen('M3B'); fidMC = fopen('M3C');
L3D_jun = []; L3A_jun = []; L3B_jun = []; L3C_jun = [];
L3D_lin = []; L3A_lin = []; L3B_lin = []; L3C_lin = [];

L3A_v2 = zeros(size(EOML3A)); 
L3B_v2 = zeros(size(EOML3B));
L3C_v2 = zeros(size(EOML3C));
L3D_v2 = zeros(size(EOML3D));
D3A_v2 = zeros(size(EOMMM23A)); 
D3B_v2 = zeros(size(EOMMM23B));
D3C_v2 = zeros(size(EOMMM23C));
D3D_v2 = zeros(size(EOMMM23D));

M3A_v2 = zeros(size(EOMMM23A)); 
M3B_v2 = zeros(size(EOMMM23B));
M3C_v2 = zeros(size(EOMMM23C));
M3D_v2 = zeros(size(EOMMM23D));

for i = 1:sys_ucc.Nocc_beta
    for j = i+1:sys_ucc.Nocc_beta
        for k = j+1:sys_ucc.Nocc_beta
            for a = 1:sys_ucc.Nvir_beta
                for b = a+1:sys_ucc.Nvir_beta
                    for c = b+1:sys_ucc.Nvir_beta

                        tlineD = fgetl(fidD); tlineA = fgetl(fidA);
                        tlineDD = fgetl(fidDD); tlineDA = fgetl(fidDA);
                        tlineMD = fgetl(fidMD); tlineMA = fgetl(fidMA);

                        L3D_jun(ct) = str2double(tlineD);
                        L3A_jun(ct) = str2double(tlineA);
                        L3D_lin(ct) = EOML3D(k,j,i,c,b,a);
                        L3A_lin(ct) = EOML3A(k,j,i,c,b,a);

                        L3A_v2(i,j,k,a,b,c) = L3A_jun(ct);
                        L3D_v2(i,j,k,a,b,c) = L3D_jun(ct); 
                        D3D_v2(a,b,c,i,j,k) = str2double(tlineDD);
                        D3A_v2(a,b,c,i,j,k) = str2double(tlineDA);
                        M3A_v2(a,b,c,i,j,k) = str2double(tlineMA);
                        M3D_v2(a,b,c,i,j,k) = str2double(tlineMD);
                        
                        if abs(L3A_jun(ct)) > 1e-6
                            fprintf('EOML3A(%d) = %4.8f\n',ct,L3A_lin(ct))
                            fprintf('L3A(%d) = %4.8f\n',ct,L3A_jun(ct))   
                        end
                        ct = ct + 1;
                    end
                end
            end
         end
     end
end

ct = 1;
for i = 1:sys_ucc.Nocc_alpha
    for j = i+1:sys_ucc.Nocc_alpha
        for k = 1:sys_ucc.Nocc_beta
            for a = 1:sys_ucc.Nvir_alpha
                for b = a+1:sys_ucc.Nvir_alpha
                    for c = 1:sys_ucc.Nvir_beta
                        tlineB = fgetl(fidB); tlineDB = fgetl(fidDB); tlineMB = fgetl(fidMB);
                        L3B_jun(ct) = str2double(tlineB);
                        L3B_lin(ct) = EOML3B(i,j,k,a,b,c);
                        L3B_v2(i,j,k,a,b,c) = L3B_jun(ct);
                        D3B_v2(a,b,c,i,j,k) = str2double(tlineDB);
                        M3B_v2(a,b,c,i,j,k) = str2double(tlineMB);
                        if abs(L3B_jun(ct)) > 1e-6
                            fprintf('EOML3B(%d) = %4.8f\n',ct,L3B_lin(ct))
                            fprintf('L3B(%d) = %4.8f\n',ct,L3B_jun(ct))   
                        end
                        ct = ct + 1;
                    end
                end
            end
         end
     end
end

ct = 1;
for i = 1:sys_ucc.Nocc_alpha
    for j = 1:sys_ucc.Nocc_beta
        for k = j+1:sys_ucc.Nocc_beta
            for a = 1:sys_ucc.Nvir_alpha
                for b = 1:sys_ucc.Nvir_beta
                    for c = b+1:sys_ucc.Nvir_beta
                        tlineC = fgetl(fidC); tlineDC = fgetl(fidDC); tlineMC = fgetl(fidMC);
                        L3C_jun(ct) = str2double(tlineC);
                        L3C_lin(ct) = EOML3C(i,j,k,a,b,c);
                        L3C_v2(i,j,k,a,b,c) = L3C_jun(ct);
                        D3C_v2(a,b,c,i,j,k) = str2double(tlineDC);
                        M3C_v2(a,b,c,i,j,k) = str2double(tlineMC);
                        if abs(L3C_jun(ct)) > 1e-6
                            fprintf('EOML3C(%d) = %4.8f\n',ct,L3C_lin(ct))
                            fprintf('L3C(%d) = %4.8f\n',ct,L3C_jun(ct))   
                        end
                        ct = ct + 1;
                    end
                end
            end
         end
     end
end
%%
iroot = 1;
[omega_crcc23A,omega_crcc23B,omega_crcc23C,omega_crcc23D] = creomucc23_debug(L3A_v2,L3B_v2,L3C_v2,L3D_v2,D3A_v2,D3B_v2,D3C_v2,D3D_v2,...
                                                                                M3A_v2,M3B_v2,M3C_v2,M3D_v2,cc_t,omega,HBar_t,sys_ucc,iroot);
%%

L3A_diff = abs(L3A_lin) - abs(L3A_jun);
[maxdiff,idx] = max(abs(L3A_diff));
norm(L3A_diff)/norm(L3A_jun)*100

L3B_diff = abs(L3B_lin) - abs(L3B_jun);
[maxdiff,idx] = max(abs(L3B_diff));
norm(L3B_diff)/norm(L3B_jun)*100

L3C_diff = abs(L3C_lin) - abs(L3C_jun);
[maxdiff,idx] = max(abs(L3C_diff));
norm(L3C_diff)/norm(L3C_jun)*100

L3D_diff = abs(L3D_lin) - abs(L3D_jun);
[maxdiff,idx] = max(abs(L3D_diff));
norm(L3D_diff)/norm(L3D_jun)*100
%%

fid = fopen('L3D');
tline = (fgetl(fid));
ct = 1;
while tline ~= -1
fprintf('L3D(%d) = %4.2f\n',ct,tline)
ct = ct + 1;
tline = (fgetl(fid));
end



%%
