clear all
clc
close all

addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab'));


%% h2o 2Re cc-pVDZ

% Integrals
addpath(genpath('/Users/karthik/Desktop/CC_matlab_tests/CAD_triples_correction/h2o-2Re/cipsi-rm-disconnected'));
load h2o-2Re-pvdz-ints.mat

Efci = -75.9516669840; % new fci number calculated by ilias
Escf = -75.587711247620774;

num_singles = 66;
num_doubles = 3349;
num_triples = 86864;
num_quadruples = 1201298;

ndet_lab = [1, 10, 100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000];
ndet = [1, 18, 167, 1399, 5664, 11350, 90880, 181579, 718316, 1390678];

nsing = [0, 2, 10, 24, 36, 38, 62, 66, 66, 66]; perc_sing = nsing/num_singles * 100;
ndoub = [0, 12, 141, 606, 1444, 1862, 2858, 3068, 3272, 3324]; perc_doub = ndoub/num_doubles * 100;
ntrip = [0, 2, 4, 298, 1550, 3164, 18750, 27542, 51812, 63590]; perc_trip = ntrip/num_triples * 100;
nquad = [0, 1, 11, 428, 2325, 5324, 44857, 81990, 256866, 383377]; perc_quad = nquad/num_quadruples * 100;

Ecipsi = [-75.587711247620774, -75.774464425666906, -75.851992681871394, -75.916671187318130, ...
          -75.937849757406752, -75.942655982203092, -75.949231100038418, -75.950248882606999, ...
          -75.951394028319214, -75.951530355513185];
      
Eextvar = Ecipsi - Escf;
      
EextprojA = [0.0000000000, -0.1867531779, -0.2642814341, -0.3289599863, ...
           -0.3501385226, -0.3549446599, -0.3615198472, -0.3625376555, ...
           -0.3636827056, -0.3638189796];
       
EextprojB = [0.0000000000, -0.1867531779, -0.2639477506, -0.3289599863, ...
            -0.3501385226, -0.3549446599, -0.3615198472, -0.3625376555, ...
            -0.3636827056, -0.3638189796];
      
EcadA = [-75.929632907722, -75.914681695303, -75.907300271504, -75.925162254484, ...
         -75.938960236188, -75.943014058700, -75.949238327047, -75.950249793421, ...
         -75.951393552679, -75.951530132110];
    
EcadB =      [-75.929632868129, -75.921405938038, -75.930703393103, -75.942028006976, ...  
              -75.947564788106, -75.948989700680, -75.950878749859, -75.951199557874, ...
              -75.951519769167, -75.951593187189];     

dat_struct.proj = 'H2O/cc-pVDZ (R/Re = 2.0)';
dat_struct.ndet_lab = ndet_lab;
dat_struct.ndet = ndet;
dat_struct.cipsi = (Ecipsi - Efci) * 1000;
dat_struct.cadA = (EcadA - Efci) * 1000;
dat_struct.cadB = (EcadB - Efci) * 1000;
dat_struct.percsing = perc_sing;
dat_struct.percdoub = perc_doub;
dat_struct.perctrip = perc_trip;
dat_struct.percquad = perc_quad;
dat_struct.fci = Efci;

workpath = '/Users/karthik/Desktop/CC_matlab_tests/CAD_triples_correction/h2o-2Re/cipsi-rm-disconnected/';

dets = [1, 10, 100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000];

EcadC = zeros(1,length(dets));

for D = 1:length(dets)
    
    run = ['ndet_',num2str(dets(D))];
    addpath(genpath([workpath,run]))
    load pspace.mat
    load matfile.mat

    nfzc = 0; nfzv = 0; 
    nact_h = 100; nact_p = 100; % BE CAREFUL ABOUT SINGLETON DIMENSIONS!
    nact_h_alpha = nact_h; nact_h_beta = nact_h; nact_p_alpha = nact_p; nact_p_beta = nact_p;
    flag_act_scheme = 0;

    sys = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,...
                               nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme);

    Ecorr_p = ucc_energy(t1a,t1b,t2a,t2b,t2c,sys);     
    Ecadcc = Ecorr_p+sys.Escf;
    fprintf('CAD-CC Energy = %4.8f\n',Ecadcc)
    fprintf('Correlation Energy = %4.8f\n',Ecadcc - sys.Escf)

    cc_t.t1a = t1a; cc_t.t1b = t1b; cc_t.t2a = t2a; cc_t.t2b = t2b; cc_t.t2c = t2c;

    [HBar_t] = build_ucc_HBar( cc_t, sys, false);

    lccopts.diis_size = 5;
    lccopts.maxit = 5000;
    lccopts.tol = 1e-8;
    lccopts.ground_shift = 0.0;

    [cc_t,lccsd_resid] = luccsd(cc_t,HBar_t,sys,lccopts);

    [EcrccP3D] = mm3_correction(cc_t,HBar_t,p_space,sys);

    E_correction = EcrccP3D - Ecadcc;

    fprintf('Moment Triples Correction = %4.12f mEh\n',E_correction * 1000);
    
    EcadC(D) = EcrccP3D;
    
    clear sys cc_t p_space t1a t1b t2a t2b t2c

end

dat_struct.cadC = (EcadC - Efci) * 1000;

%%

print_table(dat_struct,8)

%% Functions

function [EcrccP3D] = mm3_correction(cc_t,HBar_t,p_space,sys)

        %fprintf('\n')
        %fprintf('-------------------------------------------------------------------------------------\n')
        %fprintf('Performing correction for ground state -\n')
        
        Ecorr_p = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);
                
        [deltaA,deltaB,deltaC,deltaD] = crucc_P(cc_t,HBar_t,p_space,sys);
        Ecorr_crccP3A = Ecorr_p + deltaA;
        Ecorr_crccP3B = Ecorr_p + deltaB;
        Ecorr_crccP3C = Ecorr_p + deltaC;
        Ecorr_crccP3D = Ecorr_p + deltaD;
        EcrccP3A(1) = sys.Escf + Ecorr_crccP3A;
        EcrccP3B(1) = sys.Escf + Ecorr_crccP3B;
        EcrccP3C(1) = sys.Escf + Ecorr_crccP3C;
        EcrccP3D(1) = sys.Escf + Ecorr_crccP3D;
        fprintf('\n')
        fprintf('CR-UCC(P,3)_A = %4.12f Eh     Ecorr_A = %4.12f Eh     Delta_A = %4.12f Eh\n',EcrccP3A(1),Ecorr_crccP3A,deltaA)
        fprintf('CR-UCC(P,3)_B = %4.12f Eh     Ecorr_B = %4.12f Eh     Delta_B = %4.12f Eh\n',EcrccP3B(1),Ecorr_crccP3B,deltaB)
        fprintf('CR-UCC(P,3)_C = %4.12f Eh     Ecorr_C = %4.12f Eh     Delta_C = %4.12f Eh\n',EcrccP3C(1),Ecorr_crccP3C,deltaC)
        fprintf('CR-UCC(P,3)_D = %4.12f Eh     Ecorr_D = %4.12f Eh     Delta_D = %4.12f Eh\n',EcrccP3D(1),Ecorr_crccP3D,deltaD)
        fprintf('-------------------------------------------------------------------------------------\n')

end


 function [] = print_table(dat_struct, N)


    proj_name = dat_struct.proj;
    ndet_lab = dat_struct.ndet_lab; 
    ndet = dat_struct.ndet;
    cipsi = dat_struct.cipsi;
    cadA = dat_struct.cadA;
    cadB = dat_struct.cadB;
    cadC = dat_struct.cadC;
    percsing = dat_struct.percsing;
    percdoub = dat_struct.percdoub;
    perctrip = dat_struct.perctrip;
    percquad = dat_struct.percquad;
    fci = dat_struct.fci;


    fprintf('\n==============================================++%s++=====================================================\n',proj_name)
    fprintf('# Dets. Input /\n # Dets. In W.F.    %% S / %% D / %% T / %% Q        CIPSI        ec-CC (A)     ec-CC (B)    ec-CC (B;3)_TBA\n')
    for i = 1:length(ndet)
        fprintf('\n%6.0d / %6.0d    %4.2f /%4.2f /%4.2f / %4.2f    %4.6f     %4.6f     %4.6f     %4.6f',...
                 ndet_lab(i), ndet(i),...
                 round(percsing(i),3),round(percdoub(i),3),round(perctrip(i),3),round(percquad(i),3),...
                 round(cipsi(i),N),...
                 round(cadA(i),N),...
                 round(cadB(i),N),...
                 round(cadC(i),N))
    end
    fprintf('\n\nFCI total energy = %4.6f Eh\n',round(fci,8))
    fprintf('================================================================================================================\n')
 end