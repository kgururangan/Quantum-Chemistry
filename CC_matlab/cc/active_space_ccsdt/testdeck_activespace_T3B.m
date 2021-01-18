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

t1a = 100*cc_t.t1a; t1b = 100*cc_t.t1b;
t2a = 1000*cc_t.t2a; t2b = 1000*cc_t.t2b; t2c = 1000*cc_t.t2c;
t3a = 10000*cc_t.t3a; t3b = 10000*cc_t.t3b; t3c = 10000*cc_t.t3c; t3d = 10000*cc_t.t3d;

cc_t.t1a = t1a; cc_t.t1b = t1b;
cc_t.t2a = t2a; cc_t.t2b = t2b; cc_t.t2c = t2c;
cc_t.t3a = t3a; cc_t.t3b = t3b; cc_t.t3c = t3c; cc_t.t3d = t3d;

[HBar_t] = build_ucc_HBar(cc_t,sys,false);

H1A = HBar_t.H1A; H1B = HBar_t.H1B;
H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

%% Ficticious T vectors

Noa  = sys.Nocc_alpha; Nua = sys.Nvir_alpha;
Nob = sys.Nocc_beta; Nub = sys.Nvir_beta;

t1a = rand(Nua,Noa);
t1b = rand(Nub,Nob);
t2a = rand(Nua,Nua,Noa,Noa);
t2b = rand(Nua,Nub,Noa,Nob);
t2c = rand(Nub,Nub,Nob,Nob);
t3a = rand(Nua,Nua,Nua,Noa,Noa,Noa);
t3b = rand(Nua,Nua,Nub,Noa,Noa,Nob);
t3c = rand(Nua,Nub,Nub,Noa,Nob,Nob);
t3d = rand(Nub,Nub,Nub,Nob,Nob,Nob);

for a = 1:Nua
    for b = a+1:Nua
        for i = 1:Noa
            for j = i+1:Noa
                t2a(b,a,i,j) = -t2a(a,b,i,j);
                t2a(a,b,j,i) = -t2a(a,b,i,j);
                t2a(b,a,j,i) = t2a(a,b,i,j);
            end
        end
    end
end

for a = 1:Nub
    for b = a+1:Nub
        for i = 1:Nob
            for j = i+1:Nob
                t2c(b,a,i,j) = -t2c(a,b,i,j);
                t2c(a,b,j,i) = -t2c(a,b,i,j);
                t2c(b,a,j,i) = t2c(a,b,i,j);
            end
        end
    end
end

for a = 1:Nua
    for b = a+1:Nua
        for c = b+1:Nub
            for i = 1:Noa
                for j = i+1:Noa
                    for k = j+1:Noa
                            t3a(a,b,c,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(a,b,c,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(a,b,c,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(a,b,c,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(a,b,c,k,j,i) = -t3a(a,b,c,i,j,k);
                            
                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(b,a,c,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(b,a,c,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(b,a,c,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(b,a,c,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(a,c,b,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(a,c,b,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(a,c,b,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(a,c,b,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(c,b,a,i,j,k) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,k,i,j) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,j,k,i) = -t3a(a,b,c,i,j,k);
                            t3a(c,b,a,i,k,j) = t3a(a,b,c,i,j,k);
                            t3a(c,b,a,j,i,k) = t3a(a,b,c,i,j,k);
                            t3a(c,b,a,k,j,i) = t3a(a,b,c,i,j,k);
                            
                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(b,c,a,i,j,k) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(b,c,a,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(b,c,a,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(b,c,a,k,j,i) = -t3a(a,b,c,i,j,k);
                            
                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3a(c,a,b,i,j,k) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,k,i,j) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,j,k,i) = t3a(a,b,c,i,j,k);
                            t3a(c,a,b,i,k,j) = -t3a(a,b,c,i,j,k);
                            t3a(c,a,b,j,i,k) = -t3a(a,b,c,i,j,k);
                            t3a(c,a,b,k,j,i) = -t3a(a,b,c,i,j,k);
                    end
                end
            end
        end
	end
end

for a = 1:Nua
    for b = a+1:Nua
        for c = b+1:Nub
            for i = 1:Noa
                for j = i+1:Noa
                    for k = j+1:Noa
                            t3b(a,b,c,j,i,k) = -t3b(a,b,c,i,j,k);
                            t3b(b,a,c,i,j,k) = t3b(a,b,c,i,j,k);
                            t3b(b,a,c,j,i,k) = t3b(a,b,c,i,j,k);
                    end
                end
            end
        end
	end
end

for a = 1:Nua
    for b = 1:Nub
        for c = b+1:Nub
            for i = 1:Noa
                for j = 1:Nob
                    for k = j+1:Nob
                            t3c(a,b,c,i,k,j) = -t3c(a,b,c,i,j,k);
                            t3c(a,c,b,i,j,k) = -t3c(a,b,c,i,j,k);
                            t3c(a,c,b,i,k,j) = t3c(a,b,c,i,j,k);
                    end
                end
            end
        end
	end
end

for a = 1:Nub
    for b = a+1:Nub
        for c = b+1:Nub
            for i = 1:Nob
                for j = i+1:Nob
                    for k = j+1:Nob
                            t3d(a,b,c,k,i,j) = t3d(a,b,c,i,j,k);
                            t3d(a,b,c,j,k,i) = t3d(a,b,c,i,j,k);
                            t3d(a,b,c,i,k,j) = -t3d(a,b,c,i,j,k);
                            t3d(a,b,c,j,i,k) = -t3d(a,b,c,i,j,k);
                            t3d(a,b,c,k,j,i) = -t3d(a,b,c,i,j,k);
                            
                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(b,a,c,i,j,k) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,k,i,j) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,j,k,i) = -t3d(a,b,c,i,j,k);
                            t3d(b,a,c,i,k,j) = t3d(a,b,c,i,j,k);
                            t3d(b,a,c,j,i,k) = t3d(a,b,c,i,j,k);
                            t3d(b,a,c,k,j,i) = t3d(a,b,c,i,j,k);
                            
                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(a,c,b,i,j,k) = -t3d(a,b,c,i,j,k);
                            t3d(a,c,b,k,i,j) = -t3d(a,b,c,i,j,k);
                            t3d(a,c,b,j,k,i) = -t3d(a,b,c,i,j,k);
                            t3d(a,c,b,i,k,j) = t3d(a,b,c,i,j,k);
                            t3d(a,c,b,j,i,k) = t3d(a,b,c,i,j,k);
                            t3d(a,c,b,k,j,i) = t3d(a,b,c,i,j,k);
                            
                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(c,b,a,i,j,k) = -t3d(a,b,c,i,j,k);
                            t3d(c,b,a,k,i,j) = -t3d(a,b,c,i,j,k);
                            t3d(c,b,a,j,k,i) = -t3d(a,b,c,i,j,k);
                            t3d(c,b,a,i,k,j) = t3d(a,b,c,i,j,k);
                            t3d(c,b,a,j,i,k) = t3d(a,b,c,i,j,k);
                            t3d(c,b,a,k,j,i) = t3d(a,b,c,i,j,k);
                            
                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(b,c,a,i,j,k) = t3d(a,b,c,i,j,k);
                            t3d(b,c,a,k,i,j) = t3d(a,b,c,i,j,k);
                            t3d(b,c,a,j,k,i) = t3d(a,b,c,i,j,k);
                            t3d(b,c,a,i,k,j) = -t3d(a,b,c,i,j,k);
                            t3d(b,c,a,j,i,k) = -t3d(a,b,c,i,j,k);
                            t3d(b,c,a,k,j,i) = -t3d(a,b,c,i,j,k);
                            
                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3d(c,a,b,i,j,k) = t3d(a,b,c,i,j,k);
                            t3d(c,a,b,k,i,j) = t3d(a,b,c,i,j,k);
                            t3d(c,a,b,j,k,i) = t3d(a,b,c,i,j,k);
                            t3d(c,a,b,i,k,j) = -t3d(a,b,c,i,j,k);
                            t3d(c,a,b,j,i,k) = -t3d(a,b,c,i,j,k);
                            t3d(c,a,b,k,j,i) = -t3d(a,b,c,i,j,k);
                    end
                end
            end
        end
	end
end

cc_t.t1a = t1a; cc_t.t1b = t1b;
cc_t.t2a = t2a; cc_t.t2b = t2b; cc_t.t2c = t2c;
cc_t.t3a = t3a; cc_t.t3b = t3b; cc_t.t3c = t3c; cc_t.t3d = t3d; 

%%

[T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);

D3 = -einsum_kg(H1A.vv(PA,pA),T3B.PpPHHH,'Ae,BeCIJK->ABCIJK')...
+einsum_kg(H1A.vv(PA,PA),T3B.PPPHHH,'AE,EBCIJK->ABCIJK');
D3 = D3 - permute(D3,[2,1,3,4,5,6]);

D3_full = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
D3_full = D3_full - permute(D3_full,[2,1,3,4,5,6]);

get_error(D3,D3_full(PA,PA,PB,HA,HA,HB))


%% t3B - projection 1 (PPPHHH) [DONE]

clc

out_proj = 'ABCIJK';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
d1 = {'-h1A(mI),t3b(ABCmJK)'};
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
% <ijkabc | (He^T)_C | 0>
[X3B,VT3A,VT3B] = build_t3b(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,cc_t.t3a,cc_t.t3b,cc_t.t3c,cc_t.t3d,sys,0.0);
 % <IJKABC | (He^T)_C | 0>
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

out_proj = 'ABCIJk';
fun = @(x,y) hbar_term_wrap(x,out_proj,y);

% D1 = -einsum_kg(H1A.oo,t3b,'mi,abcmjk->abcijk');
d1 = {'-h1A(MI),t3b(ABCMJk)'};
fun(d1,'D1')

% D2 = -einsum_kg(H1B.oo,t3b,'mk,abcijm->abcijk');
d2 = {'-h1B(mk),t3b(ABCIJm)'};
fun(d2,'D2')

% D3 = einsum_kg(H1A.vv,t3b,'ae,ebcijk->abcijk');
d3 = {'h1A(Ae),t3b(eBCIJk)'};
fun(d3,'D3')

% D4 = einsum_kg(H1B.vv,t3b,'ce,abeijk->abcijk');
d4 = {'h1B(Ce),t3b(ABeIJk)'};
fun(d4,'D4')

% D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
d5 = {'h2A(mnIJ),t3b(ABCmnk)'};
fun(d5,'D5')

% D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
d6 = {'h2B(mnJk),t3b(ABCImn)'};
fun(d6,'D6')

% D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
d7 = {'h2A(ABef),t3b(efCIJk)'};
fun(d7,'D7')

% D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk'); 
d8 = {'h2B(BCef),t3b(AefIJk)'};
fun(d8,'D8')

% D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');  
d9 = {'h2A(AmIe),t3b(eBCmJk)'};
fun(d9,'D9')

% D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk'); 
d10 = {'h2B(AmIe),t3c(BCeJkm)'};
fun(d10,'D10')

% D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
d11 = {'h2B(mCek),t3a(ABeIJm)'};
fun(d11,'D11')

% D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
d12 = {'h2C(Cmke),t3b(ABeIJm)'};
fun(d12,'D12')

% D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
d13 = {'-h2B(Amek),t3b(eBCIJm)'};
fun(d13,'D13')

% D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
d14 = {'-h2B(mCIe),t3b(ABemJk)'};
fun(d14,'D14')

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

%%
function [] = hbar_term_wrap(arr,out_proj,label)
    write_term = @(x,y) write_einsum(x,out_proj,y);
    [D] = ccsdt_act_HBarT3_simple(arr);
    write_term(D,label)
end
