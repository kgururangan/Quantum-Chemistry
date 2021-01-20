function [cc_t,Ecorr] = uccsdt2(sys,opts,T_guess,flag_test)

    diis_size = opts.diis_size;
    maxit = opts.maxit;
    tol = opts.tol;
    shift = opts.shift;

    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

    if nargin < 3 || isempty(T_guess)
        T = zeros(sys.triples_dim,1);
    else
        T = T_guess;
    end
    
    tic_Start = tic;
    fprintf('\n==================================++Entering UCCSDT-2 Routine++=================================\n')

    % Initialize T1, T2, and T3 DIIS containers
    
    T_list = zeros(sys.triples_dim,diis_size);
    T_resid_list = zeros(sys.triples_dim,diis_size);

    szt1a = [sys.Nvir_alpha, sys.Nocc_alpha];
    szt1b = [sys.Nvir_beta, sys.Nocc_beta];
    szt2a = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nocc_alpha, sys.Nocc_alpha];
    szt2b = [sys.Nvir_alpha, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_beta];
    szt2c = [sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_beta, sys.Nocc_beta];
    szt3a = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nvir_alpha, sys.Nocc_alpha, sys.Nocc_alpha, sys.Nocc_alpha];
    szt3b = [sys.Nvir_alpha, sys.Nvir_alpha, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_alpha, sys.Nocc_beta];
    szt3c = [sys.Nvir_alpha, sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_alpha, sys.Nocc_beta, sys.Nocc_beta];
    szt3d = [sys.Nvir_beta, sys.Nvir_beta, sys.Nvir_beta, sys.Nocc_beta, sys.Nocc_beta, sys.Nocc_beta];

    t1a = reshape(T(sys.posv{1}),szt1a);
    t1b = reshape(T(sys.posv{2}),szt1b);
    t2a = reshape(T(sys.posv{3}),szt2a);
    t2b = reshape(T(sys.posv{4}),szt2b);
    t2c = reshape(T(sys.posv{5}),szt2c);
    t3a = reshape(T(sys.posv{6}),szt3a);
    t3b = reshape(T(sys.posv{7}),szt3b);
    t3c = reshape(T(sys.posv{8}),szt3c);
    t3d = reshape(T(sys.posv{9}),szt3d);

    cc_t = make_cct_struct(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d);
    [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
    
    % Jacobi/DIIS iterations
    it_micro = 0; flag_conv = 0; it_macro = 0; Ecorr_old = 0.0;
    fprintf('\nDIIS Cycle - %d',it_macro)
    while it_micro < maxit
        
        tic
        
        % store old T and get current diis dimensions
        T_old = T;
        
        t1a = reshape(T(sys.posv{1}),szt1a);
        t1b = reshape(T(sys.posv{2}),szt1b);
        t2a = reshape(T(sys.posv{3}),szt2a);
        t2b = reshape(T(sys.posv{4}),szt2b);
        t2c = reshape(T(sys.posv{5}),szt2c);
        t3a = reshape(T(sys.posv{6}),szt3a);
        %t3a = zero_t3_outside_act(t3a,2,'A',sys);
        t3b = reshape(T(sys.posv{7}),szt3b);
        %t3b = zero_t3_outside_act(t3b,2,'B',sys);
        t3c = reshape(T(sys.posv{8}),szt3c);
        %t3c = zero_t3_outside_act(t3c,2,'C',sys);
        t3d = reshape(T(sys.posv{9}),szt3d);
        %t3d = zero_t3_outside_act(t3d,2,'D',sys);

        % CC correlation energy
        Ecorr = ucc_energy(t1a,t1b,t2a,t2b,t2c,sys);
       
        % update t1 and t2 by Jacobi                          
        t1a = update_t1a_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        t1b = update_t1b_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        t2a = update_t2a_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        t2b = update_t2b_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        t2c = update_t2c_ccsdt(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
        
        % closed shell
        %t1b = t1a
        %t2c = t2a

        if ~flag_test

            [HBar_t] = build_ccsd_HBar_intermediates(t1a,t1b,t2a,t2b,t2c,sys);

            T3A.PPPHHH = update_t3a_proj1_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3A.PPpHHH = update_t3a_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3A.PPPhHH = update_t3a_proj3_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3A.PPphHH = update_t3a_proj4_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);

            T3B.PPPHHH = update_t3b_proj1_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3B.PPpHHH = update_t3b_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3B.PpPHHH = update_t3b_proj3_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3B.PPPHHh = update_t3b_proj4_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3B.PPPhHH = update_t3b_proj5_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3B.PPphHH = update_t3b_proj6_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3B.PPpHHh = update_t3b_proj7_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3B.PpPhHH = update_t3b_proj8_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3B.PpPHHh = update_t3b_proj9_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);

            T3C.PPPHHH = update_t3c_proj1_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3C.PPpHHH = update_t3c_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3C.pPPHHH = update_t3c_proj3_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3C.PPPHhH = update_t3c_proj4_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3C.PPPhHH = update_t3c_proj5_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3C.PPpHhH = update_t3c_proj6_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3C.PPphHH = update_t3c_proj7_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3C.pPPHhH = update_t3c_proj8_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3C.pPPhHH = update_t3c_proj9_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);

            T3D.PPPHHH = update_t3d_proj1_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3D.PPpHHH = update_t3d_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3D.PPPhHH = update_t3d_proj3_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            T3D.PPphHH = update_t3d_proj4_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);

            [t3a] = make_act_struct_into_t3(T3A,'A',sys);
            [t3b] = make_act_struct_into_t3(T3B,'B',sys);
            [t3c] = make_act_struct_into_t3(T3C,'C',sys);
            [t3d] = make_act_struct_into_t3(T3D,'D',sys);

        else

            [HBar_t] = build_ccsd_HBar_intermediates(t1a,t1b,t2a,t2b,t2c,sys);

            cc_t = make_cct_struct(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d);
            [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
            t3a = update_t3a(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
            t3a = zero_t3_outside_act(t3a,2,'A',sys);

            t3a_proj1 = update_t3a_proj1_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3a_proj2 = update_t3a_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3a_proj3 = update_t3a_proj3_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3a_proj4 = update_t3a_proj4_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            fprintf('\nError in t3a_proj1 = %4.15f\n',get_error(t3a_proj1,t3a(PA,PA,PA,HA,HA,HA)))
            fprintf('Error in t3a_proj2 = %4.15f\n',get_error(t3a_proj2,t3a(PA,PA,pA,HA,HA,HA)))
            fprintf('Error in t3a_proj3 = %4.15f\n',get_error(t3a_proj3,t3a(PA,PA,PA,hA,HA,HA)))
            fprintf('Error in t3a_proj4 = %4.15f\n',get_error(t3a_proj4,t3a(PA,PA,pA,hA,HA,HA))) 
  
            % adding these zero statements perfore populating with projections returns us to our original
            % CCSDt(II) energy
%            t3a_temp = t3a;
%            t3a = zeros(szt3a);
            t3a(PA,PA,PA,HA,HA,HA) = t3a_proj1;
            t3a(PA,PA,pA,HA,HA,HA) = t3a_proj2;
            t3a(PA,PA,PA,hA,HA,HA) = t3a_proj3;
            t3a(PA,PA,pA,hA,HA,HA) = t3a_proj4; 
%             t3a(PA,PA,PA,HA,HA,HA) = t3a_temp(PA,PA,PA,HA,HA,HA);
%             t3a(PA,PA,pA,HA,HA,HA) = t3a_temp(PA,PA,pA,HA,HA,HA);
%             t3a(PA,PA,PA,hA,HA,HA) = t3a_temp(PA,PA,PA,hA,HA,HA);
%             t3a(PA,PA,pA,hA,HA,HA) = t3a_temp(PA,PA,pA,hA,HA,HA);
%             get_error(t3a,t3a_temp)

            T3A.PPPHHH = t3a_proj1; T3A.PPpHHH = t3a_proj2; T3A.PPPhHH = t3a_proj3; T3A.PPphHH = t3a_proj4;
            [t3a_act] = make_act_struct_into_t3(T3A,'A',sys);
            get_error(t3a,t3a_act)

            cc_t = make_cct_struct(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d);
            [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
            t3b = update_t3b(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
            t3b = zero_t3_outside_act(t3b,2,'B',sys);

            t3b_proj1 = update_t3b_proj1_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3b_proj2 = update_t3b_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3b_proj3 = update_t3b_proj3_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3b_proj4 = update_t3b_proj4_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3b_proj5 = update_t3b_proj5_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3b_proj6 = update_t3b_proj6_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3b_proj7 = update_t3b_proj7_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3b_proj8 = update_t3b_proj8_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3b_proj9 = update_t3b_proj9_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            fprintf('\nError in t3b_proj1 = %4.15f\n',get_error(t3b_proj1,t3b(PA,PA,PB,HA,HA,HB)))
            fprintf('Error in t3b_proj2 = %4.15f\n',get_error(t3b_proj2,t3b(PA,PA,pB,HA,HA,HB)))
            fprintf('Error in t3b_proj3 = %4.15f\n',get_error(t3b_proj3,t3b(PA,pA,PB,HA,HA,HB)))
            fprintf('Error in t3b_proj4 = %4.15f\n',get_error(t3b_proj4,t3b(PA,PA,PB,HA,HA,hB)))
            fprintf('Error in t3b_proj5 = %4.15f\n',get_error(t3b_proj5,t3b(PA,PA,PB,hA,HA,HB)))
            fprintf('Error in t3b_proj6 = %4.15f\n',get_error(t3b_proj6,t3b(PA,PA,pB,hA,HA,HB)))
            fprintf('Error in t3b_proj7 = %4.15f\n',get_error(t3b_proj7,t3b(PA,PA,pB,HA,HA,hB)))
            fprintf('Error in t3b_proj8 = %4.15f\n',get_error(t3b_proj8,t3b(PA,pA,PB,hA,HA,HB)))
            fprintf('Error in t3b_proj9 = %4.15f\n',get_error(t3b_proj9,t3b(PA,pA,PB,HA,HA,hB)))

%            t3b_temp = t3b;
%            t3b = zeros(szt3b);
            t3b(PA,PA,PB,HA,HA,HB) = t3b_proj1; 
            t3b(PA,PA,pB,HA,HA,HB) = t3b_proj2; 
            t3b(PA,pA,PB,HA,HA,HB) = t3b_proj3;
            t3b(PA,PA,PB,HA,HA,hB) = t3b_proj4; 
            t3b(PA,PA,PB,hA,HA,HB) = t3b_proj5; 
            t3b(PA,PA,pB,hA,HA,HB) = t3b_proj6;
            t3b(PA,PA,pB,HA,HA,hB) = t3b_proj7; 
            t3b(PA,pA,PB,hA,HA,HB) = t3b_proj8; 
            t3b(PA,pA,PB,HA,HA,hB) = t3b_proj9;
%             t3b(PA,PA,PB,HA,HA,HB) = t3b_temp(PA,PA,PB,HA,HA,HB); 
%             t3b(PA,PA,pB,HA,HA,HB) = t3b_temp(PA,PA,pB,HA,HA,HB); 
%             t3b(PA,pA,PB,HA,HA,HB) = t3b_temp(PA,pA,PB,HA,HA,HB);
%             t3b(PA,PA,PB,HA,HA,hB) = t3b_temp(PA,PA,PB,HA,HA,hB); 
%             t3b(PA,PA,PB,hA,HA,HB) = t3b_temp(PA,PA,PB,hA,HA,HB); 
%             t3b(PA,PA,pB,hA,HA,HB) = t3b_temp(PA,PA,pB,hA,HA,HB);
%             t3b(PA,PA,pB,HA,HA,hB) = t3b_temp(PA,PA,pB,HA,HA,hB); 
%             t3b(PA,pA,PB,hA,HA,HB) = t3b_temp(PA,pA,PB,hA,HA,HB); 
%             t3b(PA,pA,PB,HA,HA,hB) = t3b_temp(PA,pA,PB,HA,HA,hB);
%             get_error(t3b,t3b_temp)

            T3B.PPPHHH = t3b_proj1; T3B.PPpHHH = t3b_proj2; T3B.PpPHHH = t3b_proj3;
            T3B.PPPHHh = t3b_proj4; T3B.PPPhHH = t3b_proj5; T3B.PPphHH = t3b_proj6;
            T3B.PPpHHh = t3b_proj7; T3B.PpPhHH = t3b_proj8; T3B.PpPHHh = t3b_proj9;
            [t3b_act] = make_act_struct_into_t3(T3B,'B',sys);
            get_error(t3b,t3b_act)

            cc_t = make_cct_struct(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d);
            [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
            t3c = update_t3c(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
            t3c = zero_t3_outside_act(t3c,2,'C',sys);

            t3c_proj1 = update_t3c_proj1_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3c_proj2 = update_t3c_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3c_proj3 = update_t3c_proj3_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3c_proj4 = update_t3c_proj4_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3c_proj5 = update_t3c_proj5_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3c_proj6 = update_t3c_proj6_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3c_proj7 = update_t3c_proj7_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3c_proj8 = update_t3c_proj8_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3c_proj9 = update_t3c_proj9_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            fprintf('\nError in t3c_proj1 = %4.15f\n',get_error(t3c_proj1,t3c(PA,PB,PB,HA,HB,HB)))
            fprintf('Error in t3c_proj2 = %4.15f\n',get_error(t3c_proj2,t3c(PA,PB,pB,HA,HB,HB)))
            fprintf('Error in t3c_proj3 = %4.15f\n',get_error(t3c_proj3,t3c(pA,PB,PB,HA,HB,HB)))
            fprintf('Error in t3c_proj4 = %4.15f\n',get_error(t3c_proj4,t3c(PA,PB,PB,HA,hB,HB)))
            fprintf('Error in t3c_proj5 = %4.15f\n',get_error(t3c_proj5,t3c(PA,PB,PB,hA,HB,HB)))
            fprintf('Error in t3c_proj6 = %4.15f\n',get_error(t3c_proj6,t3c(PA,PB,pB,HA,hB,HB)))
            fprintf('Error in t3c_proj7 = %4.15f\n',get_error(t3c_proj7,t3c(PA,PB,pB,hA,HB,HB)))
            fprintf('Error in t3c_proj8 = %4.15f\n',get_error(t3c_proj8,t3c(pA,PB,PB,HA,hB,HB)))
            fprintf('Error in t3c_proj9 = %4.15f\n',get_error(t3c_proj9,t3c(pA,PB,PB,hA,HB,HB)))

%            t3c_temp = t3c;
%            t3c = zeros(szt3c);
            t3c(PA,PB,PB,HA,HB,HB) = t3c_proj1; 
            t3c(PA,PB,pB,HA,HB,HB) = t3c_proj2; 
            t3c(pA,PB,PB,HA,HB,HB) = t3c_proj3;
            t3c(PA,PB,PB,HA,hB,HB) = t3c_proj4; 
            t3c(PA,PB,PB,hA,HB,HB) = t3c_proj5; 
            t3c(PA,PB,pB,HA,hB,HB) = t3c_proj6;
            t3c(PA,PB,pB,hA,HB,HB) = t3c_proj7; 
            t3c(pA,PB,PB,HA,hB,HB) = t3c_proj8; 
            t3c(pA,PB,PB,hA,HB,HB) = t3c_proj9;
%             t3c(PA,PB,PB,HA,HB,HB) = t3c_temp(PA,PB,PB,HA,HB,HB); 
%             t3c(PA,PB,pB,HA,HB,HB) = t3c_temp(PA,PB,pB,HA,HB,HB); 
%             t3c(pA,PB,PB,HA,HB,HB) = t3c_temp(pA,PB,PB,HA,HB,HB);
%             t3c(PA,PB,PB,HA,hB,HB) = t3c_temp(PA,PB,PB,HA,hB,HB); 
%             t3c(PA,PB,PB,hA,HB,HB) = t3c_temp(PA,PB,PB,hA,HB,HB); 
%             t3c(PA,PB,pB,HA,hB,HB) = t3c_temp(PA,PB,pB,HA,hB,HB);
%             t3c(PA,PB,pB,hA,HB,HB) = t3c_temp(PA,PB,pB,hA,HB,HB); 
%             t3c(pA,PB,PB,HA,hB,HB) = t3c_temp(pA,PB,PB,HA,hB,HB); 
%             t3c(pA,PB,PB,hA,HB,HB) = t3c_temp(pA,PB,PB,hA,HB,HB);
%             get_error(t3c,t3c_temp)

            T3C.PPPHHH = t3c_proj1; T3C.PPpHHH = t3c_proj2; T3C.pPPHHH = t3c_proj3; 
            T3C.PPPHhH = t3c_proj4; T3C.PPPhHH = t3c_proj5; T3C.PPpHhH = t3c_proj6;
            T3C.PPphHH = t3c_proj7; T3C.pPPHhH = t3c_proj8; T3C.pPPhHH = t3c_proj9;
            [t3c_act] = make_act_struct_into_t3(T3C,'C',sys);
            get_error(t3c,t3c_act)

            cc_t = make_cct_struct(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d);
            [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys);
            t3d = update_t3d(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys,shift);
            t3d = zero_t3_outside_act(t3d,2,'D',sys);

            t3d_proj1 = update_t3d_proj1_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3d_proj2 = update_t3d_proj2_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3d_proj3 = update_t3d_proj3_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            t3d_proj4 = update_t3d_proj4_ccsdt2_v2(t1a, t1b, t2a, t2b, t2c, T3A, T3B, T3C, T3D, HBar_t, sys, shift);
            fprintf('\nError in t3d_proj1 = %4.15f\n',get_error(t3d_proj1,t3d(PB,PB,PB,HB,HB,HB)))
            fprintf('Error in t3d_proj2 = %4.15f\n',get_error(t3d_proj2,t3d(PB,PB,pB,HB,HB,HB)))
            fprintf('Error in t3d_proj3 = %4.15f\n',get_error(t3d_proj3,t3d(PB,PB,PB,hB,HB,HB)))
            fprintf('Error in t3d_proj4 = %4.15f\n',get_error(t3d_proj4,t3d(PB,PB,pB,hB,HB,HB))) 

%            t3d_temp = t3d;
%            t3d = zeros(szt3d);
            t3d(PB,PB,PB,HB,HB,HB) = t3d_proj1;
            t3d(PB,PB,pB,HB,HB,HB) = t3d_proj2;
            t3d(PB,PB,PB,hB,HB,HB) = t3d_proj3;
            t3d(PB,PB,pB,hB,HB,HB) = t3d_proj4; 
%             t3d(PB,PB,PB,HB,HB,HB) = t3d_temp(PB,PB,PB,HB,HB,HB);
%             t3d(PB,PB,pB,HB,HB,HB) = t3d_temp(PB,PB,pB,HB,HB,HB);
%             t3d(PB,PB,PB,hB,HB,HB) = t3d_temp(PB,PB,PB,hB,HB,HB);
%             t3d(PB,PB,pB,hB,HB,HB) = t3d_temp(PB,PB,pB,hB,HB,HB); 
%             get_error(t3d,t3d_temp)

             T3D.PPPHHH = t3d_proj1; T3D.PPpHHH = t3d_proj2; T3D.PPPhHH = t3d_proj3; T3D.PPphHH = t3d_proj4;
             [t3d_act] = make_act_struct_into_t3(T3D,'D',sys);
             get_error(t3d,t3d_act)
 
%             [t3a_act] = make_act_struct_into_t3(T3A,'A',sys);
%             [t3b_act] = make_act_struct_into_t3(T3B,'B',sys);
%             [t3c_act] = make_act_struct_into_t3(T3C,'C',sys);
%             [t3d_act] = make_act_struct_into_t3(T3D,'D',sys);
% 
%             fprintf('\nError in t3d_proj1 = %4.15f\n',get_error(t3d_proj1,t3d(PB,PB,PB,HB,HB,HB)))
%             fprintf('Error in t3d_proj2 = %4.15f\n',get_error(t3d_proj2,t3d(PB,PB,pB,HB,HB,HB)))
%             fprintf('Error in t3d_proj3 = %4.15f\n',get_error(t3d_proj3,t3d(PB,PB,PB,hB,HB,HB)))
%             fprintf('Error in t3d_proj4 = %4.15f\n',get_error(t3d_proj4,t3d(PB,PB,pB,hB,HB,HB)))          
        end
     
        % store vectorized results
        T(sys.posv{1}) = t1a(:); T(sys.posv{2}) = t1b(:);
        T(sys.posv{3}) = t2a(:); T(sys.posv{4}) = t2b(:); T(sys.posv{5}) = t2c(:);
        T(sys.posv{6}) = t3a(:); T(sys.posv{7}) = t3b(:); T(sys.posv{8}) = t3c(:); T(sys.posv{9}) = t3d(:);

        % build DIIS residual
        T_resid = T - T_old;  
        
        % change in Ecorr
        deltaE = Ecorr - Ecorr_old;

        % check for exit condition
        ccsdt_resid = sqrt(mean(T_resid.^2));
        if ccsdt_resid < tol && abs(deltaE) < tol
            flag_conv = 1;
            break;
        end

        % append trial and residual vectors to lists
        T_list(:,mod(it_micro,diis_size)+1) = T;
        T_resid_list(:,mod(it_micro,diis_size)+1) = T_resid;
         
        % diis extrapolate
        if mod(it_micro,diis_size) == 0 && it_micro > 1
           it_macro = it_macro + 1;
           fprintf('\nDIIS Cycle - %d',it_macro)
           T = diis_xtrap(T_list,T_resid_list);
        end        

        % extract time per iteration in minutes and seconds
        toc_S = toc; toc_M = floor(toc_S/60); toc_S = toc_S - toc_M*60;

        fprintf('\n   Iter-%d    Residuum = %4.10f    dE = %4.10f    Ecorr = %4.10f   (%dm %.1fs)',it_micro,ccsdt_resid,deltaE,Ecorr,toc_M,toc_S);
        
        it_micro = it_micro + 1;
        Ecorr_old = Ecorr;
         
    end

    % return final amplitudes and correlation energy
    cc_t.t1a = reshape(T(sys.posv{1}),szt1a);
    cc_t.t1b = reshape(T(sys.posv{2}),szt1b);
    cc_t.t2a = reshape(T(sys.posv{3}),szt2a);
    cc_t.t2b = reshape(T(sys.posv{4}),szt2b);
    cc_t.t2c = reshape(T(sys.posv{5}),szt2c);
    cc_t.t3a = reshape(T(sys.posv{6}),szt3a);
    cc_t.t3b = reshape(T(sys.posv{7}),szt3b);
    cc_t.t3c = reshape(T(sys.posv{8}),szt3c);
    cc_t.t3d = reshape(T(sys.posv{9}),szt3d);
        
    Ecorr = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);
    
    if flag_conv == 1
        fprintf('\nUCCSDT-2 successfully converged in %d iterations (%4.2f seconds)\n',it_micro,toc(tic_Start));
        fprintf('Total Energy = %4.12f Eh     Ecorr = %4.12f Eh\n',Ecorr+sys.Escf,Ecorr);
    else
        fprintf('\nUCCSDT-2 failed to converged in %d iterations\n',maxit)
    end   
        

end

function [cc_t] = make_cct_struct(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d)
    cc_t.t1a = t1a; cc_t.t1b = t1b;
    cc_t.t2a = t2a; cc_t.t2b = t2b; cc_t.t2c = t2c;
    cc_t.t3a = t3a; cc_t.t3b = t3b; cc_t.t3c = t3c; cc_t.t3d = t3d;
end
