clear all
clc
close all

load h2o-631g-stretched

nfzc = 0; nfzv = 0; 
nact_h = 5; nact_p = 10; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
sys = build_system_cc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,nfzv,nact_h,nact_p);

maxit = 5;

Wavefunction.Dets(1,:) = [1:sys.Nelec];
Wavefunction.Energy = sys.Escf;
Wavefunction.Coef = [1];
[V_ext] = get_external_space(Wavefunction,sys);

%for it = 2:maxit
    
    [eMP2_array] = cipsi_MP2(Wavefunction,V_ext,sys);
    EMP2 = sum(eMP2_array);
    
%end


