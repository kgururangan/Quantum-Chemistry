function [HBar_conv, HBar, t1_conv, t2_conv] = check_HBar(cc_t,sys_ucc,sys_cc)

      flag_build_3body = true;

%     clc
%     close all
%     clear all
% 
%     load h2o-631g-stretched
%     nfzc = 0; nfzv = 0; nact_h = 200; nact_p = 200;
% 
%     sys_ucc = build_system_ucc(e1int,e2int,Vnuc,Nocc_a,Nocc_b);
%     sys_cc = build_system_cc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,nfzv,nact_h,nact_p);
% 
%     % spinorbital indices for alpha and beta orbitals
%     occ = 1:Nelec; unocc = [Nelec+1:2*sys_ucc.Norb] - Nelec;
%     oa = occ(1:2:end); ob = occ(2:2:end);
%     ua = unocc(1:2:end); ub = unocc(2:2:end);
% 
%     %% UCCSD
% 
%     ccopts.diis_size = 5;
%     ccopts.maxit = 100;
%     ccopts.tol = 1e-9;
%     [cc_t,Ecorr_ucc] = uccsd(sys_ucc,ccopts);
%     [t1_conv] = convert_spinint_to_spinorb({cc_t.t1a,cc_t.t1b},sys_ucc);
%     [t2_conv] = convert_spinint_to_spinorb({cc_t.t2a,cc_t.t2b,cc_t.t2c},sys_ucc);

    % Build HBars
    
    [t1_conv] = convert_spinint_to_spinorb({cc_t.t1a,cc_t.t1b},sys_ucc);
    [t2_conv] = convert_spinint_to_spinorb({cc_t.t2a,cc_t.t2b,cc_t.t2c},sys_ucc);

    [HBar] = build_HBar_debug(t1_conv,t2_conv,sys_cc);

    [HBar_t] = build_ucc_HBar( cc_t, sys_ucc, flag_build_3body );
    H1A = HBar_t.H1A; H1B = HBar_t.H1B; 
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3C = HBar_t.H3C; H3D = HBar_t.H3D;

    HBar_conv = convert_HBar_to_spinorb(HBar_t);

    %
    clc

    fprintf('\n======================================== ++ 1-Body HBar ++ ========================================\n')

    % H(mi)
    fprintf('\n===================== H(oo) =====================\n')
    err_oo = HBar_conv{1}{1,1} - HBar{1}{1,1};
    fprintf('Error in H(mi) = %4.10f\n',sum(abs(err_oo(:))))
    fprintf('=================================================\n')

    % H(ae)
    fprintf('\n===================== H(vv) =====================\n')
    err_vv = HBar_conv{1}{2,2} - HBar{1}{2,2};
    fprintf('Error in H(ae) = %4.10f\n',sum(abs(err_vv(:))))
    fprintf('=================================================\n')

    % H(me)
    fprintf('\n===================== H(mv) =====================\n')
    err_ov = HBar_conv{1}{1,2} - HBar{1}{1,2};
    fprintf('Error in H(me) = %4.10f\n',sum(abs(err_ov(:))))
    fprintf('=================================================\n')

    fprintf('\n======================================== ++ 2-Body HBar ++ ========================================\n')

    % H(amij)
    fprintf('\n==================== H(vooo) ====================\n')
    err_vooo = HBar_conv{2}{2,1,1,1} - HBar{2}{2,1,1,1};
    err_ovoo = HBar_conv{2}{1,2,1,1} - HBar{2}{1,2,1,1};
    fprintf('Error in H(vooo) = %4.10f\n',sum(abs(err_vooo(:))))
    fprintf('Error in H(ovoo) = %4.10f\n',sum(abs(err_ovoo(:))))
    fprintf('=================================================\n')

    % H(abie)
    fprintf('\n==================== H(vvov) ====================\n')
    err_vvov = HBar_conv{2}{2,2,1,2} - HBar{2}{2,2,1,2};
    err_vvvo = HBar_conv{2}{2,2,2,1} - HBar{2}{2,2,2,1};
    fprintf('Error in H(vvov) = %4.10f\n',sum(abs(err_vvov(:))))
    fprintf('Error in H(vvvo) = %4.10f\n',sum(abs(err_vvvo(:))))
    fprintf('=================================================\n')

    % H(amie)
    fprintf('\n==================== H(voov) ====================\n')
    err_voov = HBar_conv{2}{2,1,1,2} - HBar{2}{2,1,1,2};
    err_ovov = HBar_conv{2}{1,2,1,2} - HBar{2}{1,2,1,2};
    err_vovo = HBar_conv{2}{2,1,2,1} - HBar{2}{2,1,2,1};
    err_ovvo = HBar_conv{2}{1,2,2,1} - HBar{2}{1,2,2,1};
    fprintf('Error in H(voov) = %4.10f\n',sum(abs(err_voov(:))))
    fprintf('Error in H(ovov) = %4.10f\n',sum(abs(err_ovov(:))))
    fprintf('Error in H(vovo) = %4.10f\n',sum(abs(err_vovo(:))))
    fprintf('Error in H(ovvo) = %4.10f\n',sum(abs(err_ovvo(:))))
    fprintf('=================================================\n')

    % H(abef)
    fprintf('\n==================== H(vvvv) ====================\n')
    err_vvvv = HBar_conv{2}{2,2,2,2} - HBar{2}{2,2,2,2};
    fprintf('Error in H(vvvv) = %4.10f\n',sum(abs(err_vvvv(:))))
    fprintf('=================================================\n')

    % H(mnij)
    fprintf('\n==================== H(oooo) ====================\n')
    err_oooo = HBar_conv{2}{1,1,1,1} - HBar{2}{1,1,1,1};
    fprintf('Error in H(oooo) = %4.10f\n',sum(abs(err_oooo(:))))
    fprintf('=================================================\n')

    % H(amef)
    fprintf('\n==================== H(vovv) ====================\n')
    err_vovv = HBar_conv{2}{2,1,2,2} - HBar{2}{2,1,2,2};
    err_ovvv = HBar_conv{2}{1,2,2,2} - HBar{2}{1,2,2,2};
    fprintf('Error in H(vovv) = %4.10f\n',sum(abs(err_vovv(:))))
    fprintf('Error in H(ovvv) = %4.10f\n',sum(abs(err_ovvv(:))))
    fprintf('=================================================\n')

    % H(mnie)
    fprintf('\n==================== H(ooov) ====================\n')
    err_ooov = HBar_conv{2}{1,1,1,2} - HBar{2}{1,1,1,2};
    err_oovo = HBar_conv{2}{1,1,2,1} - HBar{2}{1,1,2,1};
    fprintf('Error in H(ooov) = %4.10f\n',sum(abs(err_ooov(:))))
    fprintf('Error in H(oovo) = %4.10f\n',sum(abs(err_oovo(:))))
    fprintf('=================================================\n')

    % H(mnef)
    fprintf('\n==================== H(oovv) ====================\n')
    err_oovv = HBar_conv{2}{1,1,2,2} - HBar{2}{1,1,2,2};
    fprintf('Error in H(oovv) = %4.10f\n',sum(abs(err_oovv(:))))
    fprintf('=================================================\n')
    
    fprintf('\n======================================== ++ 3-Body HBar ++ ========================================\n')
    
    % H(abmije)
    fprintf('\n=================== H(vvooov) ===================\n')
    err_vvooov = HBar_conv{3}{2,2,1,1,1,2} - HBar{3}{2,2,1,1,1,2};
    fprintf('Error in H(vvooov) = %4.10f\n',sum(abs(err_vvooov(:))))
    fprintf('=================================================\n')
    
    % H(abmief)
    fprintf('\n=================== H(vvoovv) ===================\n')
    err_vvoovv = HBar_conv{3}{2,2,1,1,2,2} - HBar{3}{2,2,1,1,2,2};
    fprintf('Error in H(vvoovv) = %4.10f\n',sum(abs(err_vvoovv(:))))
    fprintf('=================================================\n')
end

    % %% Checking the individual HBar spin-integrated parts

    % clc
    % 
    % fprintf('\n======================================== ++ 1-Body HBar ++ ========================================\n')
    % 
    % % H(mi)
    % fprintf('\n===================== H(oo) =====================\n')
    % 
    % % exact H(mi) in spinorbital form
    % Hex_oo = HBar{1}{1,1};
    % 
    % % spinorbital container
    % Hchk_oo = zeros(size(Hex_oo));
    % 
    % Hchk_oo(oa,oa) = H1A.oo;
    % Hchk_oo(ob,ob) = H1B.oo;
    % 
    % errA_oo = abs(Hchk_oo(oa,oa) - Hex_oo(oa,oa));
    % errB_oo = abs(Hchk_oo(ob,ob) - Hex_oo(ob,ob));
    % 
    % fprintf('Error in H1A(mi) = %4.10f\n',sum(errA_oo(:)))
    % fprintf('Error in H1B(mi) = %4.10f\n',sum(errB_oo(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % % H(ae)
    % fprintf('\n===================== H(vv) =====================\n')
    % 
    % % exact H(ae) in spinorbital form
    % Hex_vv = HBar{1}{2,2};
    % 
    % % spinorbital container
    % Hchk_vv = zeros(size(Hex_vv));
    % 
    % Hchk_vv(ua,ua) = H1A.vv;
    % Hchk_vv(ub,ub) = H1B.vv;
    % 
    % errA_vv = abs(Hchk_vv(ua,ua) - Hex_vv(ua,ua));
    % errB_vv = abs(Hchk_vv(ub,ub) - Hex_vv(ub,ub));
    % 
    % fprintf('Error in H1A(ae) = %4.10f\n',sum(errA_vv(:)))
    % fprintf('Error in H1B(ae) = %4.10f\n',sum(errB_vv(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % % H(me)
    % fprintf('\n===================== H(ov) =====================\n')
    % 
    % % exact H(ae) in spinorbital form
    % Hex_ov = HBar{1}{1,2};
    % 
    % % spinorbital container
    % Hchk_ov = zeros(size(Hex_ov));
    % 
    % Hchk_ov(oa,ua) = H1A.ov;
    % Hchk_ov(ob,ub) = H1B.ov;
    % 
    % errA_ov = abs(Hchk_ov(oa,ua) - Hex_ov(oa,ua));
    % errB_ov = abs(Hchk_ov(ob,ub) - Hex_ov(ob,ub));
    % 
    % fprintf('Error in H1A(me) = %4.10f\n',sum(errA_ov(:)))
    % fprintf('Error in H1B(me) = %4.10f\n',sum(errB_ov(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % fprintf('\n======================================== ++ 2-Body HBar ++ ========================================\n')
    % 
    % % H(amij)
    % fprintf('\n==================== H(vooo) ====================\n')
    % 
    % % exact H(amij) in spinorbital form
    % Hex_vooo = HBar{2}{2,1,1,1};
    % Hex_ovoo = HBar{2}{1,2,1,1};
    % 
    % % spinorbital container
    % Hchk_vooo = zeros(size(Hex_vooo));
    % Hchk_ovoo = zeros(size(Hex_ovoo));
    % 
    % % need to find a way to populate all the following cases:
    % % h(ua,oa,oa,oa), h(ub,ob,ob,ob), h(ua,ob,oa,ob), h(ub,oa,ob,oa),
    % % h(ub,oa,oa,ob), h(ua,ob,ob,oa)
    % % using spin-integrated H.vooo and H.ovoo
    % % E.g.
    % % If you need to get Hvooo(ub,oa,ob,oa), you need to permute Hovoo since
    % % Hvooo can't produce ub. The permutation is permute(H2B.ovoo,[2,1,4,3]) which is 2 swaps
    % % so a sign of +1. On the other hand, producing 
    % 
    % Hchk_vooo(ua,oa,oa,oa) = H2A.vooo;
    % Hchk_ovoo(oa,ua,oa,oa) = H2A.ovoo;
    % 
    % Hchk_vooo(ub,ob,ob,ob) = H2C.vooo;
    % Hchk_ovoo(ob,ub,ob,ob) = H2C.ovoo;
    % 
    % Hchk_vooo(ua,ob,oa,ob) = H2B.vooo;
    % Hchk_vooo(ub,oa,ob,oa) = permute(H2B.ovoo,[2,1,4,3]);
    % Hchk_vooo(ub,oa,oa,ob) = -permute(H2B.ovoo,[2,1,3,4]);
    % Hchk_vooo(ua,ob,ob,oa) = -permute(H2B.vooo,[1,2,4,3]);
    % 
    % 
    % Hchk_ovoo(oa,ub,oa,ob) = H2B.ovoo;
    % 
    % 
    % errA_vooo = abs(Hchk_vooo(ua,oa,oa,oa) - Hex_vooo(ua,oa,oa,oa));
    % errC_vooo = abs(Hchk_vooo(ub,ob,ob,ob) - Hex_vooo(ub,ob,ob,ob));
    % errB_vooo = abs(Hchk_vooo(ua,ob,oa,ob) - Hex_vooo(ua,ob,oa,ob));
    % err_vooo = abs(Hchk_vooo - Hex_vooo);
    % %err_vooo = abs(Hchk_vooo(ub,oa,oa,ob) - Hex_vooo(ub,oa,oa,ob));
    % 
    % errA_ovoo = abs(Hchk_ovoo(oa,ua,oa,oa) - Hex_ovoo(oa,ua,oa,oa));
    % errC_ovoo = abs(Hchk_ovoo(ob,ub,ob,ob) - Hex_ovoo(ob,ub,ob,ob));
    % errB_ovoo = abs(Hchk_ovoo(oa,ub,oa,ob) - Hex_ovoo(oa,ub,oa,ob));
    % 
    % fprintf('Error in H2A(amij) = %4.10f\n',sum(errA_vooo(:)))
    % fprintf('Error in H2C(amij) = %4.10f\n',sum(errC_vooo(:)))
    % fprintf('Error in H2B(amij) = %4.10f\n',sum(errB_vooo(:)))
    % fprintf('Total error in H(amij) = %4.10f\n',sum(err_vooo(:)))
    % fprintf('\n')
    % fprintf('Error in H2A(maij) = %4.10f\n',sum(errA_ovoo(:)))
    % fprintf('Error in H2C(maij) = %4.10f\n',sum(errC_ovoo(:)))
    % fprintf('Error in H2B(maji) = %4.10f\n',sum(errB_ovoo(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % % H(abie)
    % fprintf('\n==================== H(vvov) ====================\n')
    % 
    % % exact H(abie) in spinorbital form
    % Hex_vvov = HBar{2}{2,2,1,2};
    % Hex_vvvo = HBar{2}{2,2,2,1};
    % 
    % % spinorbital container
    % Hchk_vvov = zeros(size(Hex_vvov));
    % Hchk_vvvo = zeros(size(Hex_vvvo));
    % 
    % Hchk_vvov(ua,ua,oa,ua) = H2A.vvov;
    % Hchk_vvvo(ua,ua,ua,oa) = H2A.vvvo;
    % 
    % Hchk_vvov(ub,ub,ob,ub) = H2C.vvov;
    % Hchk_vvvo(ub,ub,ub,ob) = H2C.vvvo;
    % 
    % Hchk_vvov(ua,ub,oa,ub) = H2B.vvov;
    % Hchk_vvvo(ua,ub,ua,ob) = H2B.vvvo;
    % 
    % errA_vvov = abs(Hchk_vvov(ua,ua,oa,ua) - Hex_vvov(ua,ua,oa,ua));
    % errC_vvov = abs(Hchk_vvov(ub,ub,ob,ub) - Hex_vvov(ub,ub,ob,ub));
    % errB_vvov = abs(Hchk_vvov(ua,ub,oa,ub) - Hex_vvov(ua,ub,oa,ub));
    % 
    % errA_vvvo = abs(Hchk_vvvo(ua,ua,ua,oa) - Hex_vvvo(ua,ua,ua,oa));
    % errC_vvvo = abs(Hchk_vvvo(ub,ub,ub,ob) - Hex_vvvo(ub,ub,ub,ob));
    % errB_vvvo = abs(Hchk_vvvo(ua,ub,ua,ob) - Hex_vvvo(ua,ub,ua,ob));
    % 
    % fprintf('Error in H2A(abie) = %4.10f\n',sum(errA_vvov(:)))
    % fprintf('Error in H2C(abie) = %4.10f\n',sum(errC_vvov(:)))
    % fprintf('Error in H2B(abie) = %4.10f\n',sum(errB_vvov(:)))
    % fprintf('\n')
    % fprintf('Error in H2A(abei) = %4.10f\n',sum(errA_vvvo(:)))
    % fprintf('Error in H2C(abei) = %4.10f\n',sum(errC_vvvo(:)))
    % fprintf('Error in H2B(abei) = %4.10f\n',sum(errB_vvvo(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % 
    % % H(mnij)
    % fprintf('\n==================== H(oooo) ====================\n')
    % 
    % % exact H(mnij) in spinorbital form
    % Hex_oooo = HBar{2}{1,1,1,1};
    % 
    % % spinorbital container
    % Hchk_oooo = zeros(size(Hex_oooo));
    % 
    % Hchk_oooo(oa,oa,oa,oa) = H2A.oooo;
    % Hchk_oooo(oa,ob,oa,ob) = H2B.oooo;
    % Hchk_oooo(ob,oa,oa,ob) = -permute(H2B.oooo,[2,1,3,4]);
    % Hchk_oooo(oa,ob,ob,oa) = -permute(H2B.oooo,[1,2,4,3]);
    % Hchk_oooo(ob,oa,ob,oa) = permute(H2B.oooo,[2,1,4,3]);
    % Hchk_oooo(ob,ob,ob,ob) = H2C.oooo;
    % 
    % errA_oooo = abs(Hchk_oooo(oa,oa,oa,oa) - Hex_oooo(oa,oa,oa,oa));
    % errB_oooo = abs(Hchk_oooo(oa,ob,oa,ob) - Hex_oooo(oa,ob,oa,ob));
    % errC_oooo = abs(Hchk_oooo(ob,ob,ob,ob) - Hex_oooo(ob,ob,ob,ob));
    % err_oooo = abs(Hchk_oooo - Hex_oooo);
    % 
    % fprintf('Error in H2A(mnij) = %4.10f\n',sum(errA_oooo(:)))
    % fprintf('Error in H2C(mnij) = %4.10f\n',sum(errC_oooo(:)))
    % fprintf('Error in H2B(mnij) = %4.10f\n',sum(errB_oooo(:)))
    % fprintf('Total error in H(mnij) = %4.10f\n',sum(err_oooo(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % 
    % % H(abef)
    % fprintf('\n==================== H(vvvv) ====================\n')
    % 
    % % exact H(abef) in spinorbital form
    % Hex_vvvv = HBar{2}{2,2,2,2};
    % 
    % % spinorbital container
    % Hchk_vvvv = zeros(size(Hex_vvvv));
    % 
    % Hchk_vvvv(ua,ua,ua,ua) = H2A.vvvv;
    % Hchk_vvvv(ua,ub,ua,ub) = H2B.vvvv;
    % Hchk_vvvv(ub,ub,ub,ub) = H2C.vvvv;
    % 
    % errA_vvvv = abs(Hchk_vvvv(ua,ua,ua,ua) - Hex_vvvv(ua,ua,ua,ua));
    % errB_vvvv = abs(Hchk_vvvv(ua,ub,ua,ub) - Hex_vvvv(ua,ub,ua,ub));
    % errC_vvvv = abs(Hchk_vvvv(ub,ub,ub,ub) - Hex_vvvv(ub,ub,ub,ub));
    % 
    % fprintf('Error in H2A(abef) = %4.10f\n',sum(errA_vvvv(:)))
    % fprintf('Error in H2C(abef) = %4.10f\n',sum(errC_vvvv(:)))
    % fprintf('Error in H2B(abef) = %4.10f\n',sum(errB_vvvv(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % % H(amie)
    % fprintf('\n==================== H(voov) ====================\n')
    % 
    % % exact H(amie) in spinorbital form
    % Hex_voov = HBar{2}{2,1,1,2};
    % Hex_ovov = HBar{2}{1,2,1,2};
    % Hex_vovo = HBar{2}{2,1,2,1};
    % Hex_ovvo = HBar{2}{1,2,2,1};
    % 
    % % spinorbital container
    % Hchk_voov = zeros(size(Hex_voov));
    % Hchk_ovov = zeros(size(Hex_ovov));
    % Hchk_vovo = zeros(size(Hex_vovo));
    % Hchk_ovvo = zeros(size(Hex_ovvo));
    % 
    % Hchk_voov(ua,oa,oa,ua) = H2A.voov;
    % Hchk_voov(ua,ob,oa,ub) = H2B.voov;
    % Hchk_ovov(oa,ub,oa,ub) = H2B.ovov;
    % Hchk_vovo(ua,ob,ua,ob) = H2B.vovo;
    % Hchk_ovvo(oa,ub,ua,ob) = H2B.ovvo;
    % Hchk_voov(ub,ob,ob,ub) = H2C.voov;
    % 
    % errA_voov = abs(Hchk_voov(ua,oa,oa,ua) - Hex_voov(ua,oa,oa,ua));
    % errB_voov = abs(Hchk_voov(ua,ob,oa,ub) - Hex_voov(ua,ob,oa,ub));
    % errB_ovov = abs(Hchk_ovov(oa,ub,oa,ub) - Hex_ovov(oa,ub,oa,ub));
    % errB_vovo = abs(Hchk_vovo(ua,ob,ua,ob) - Hex_vovo(ua,ob,ua,ob));
    % errB_ovvo = abs(Hchk_ovvo(oa,ub,ua,ob) - Hex_ovvo(oa,ub,ua,ob));
    % errC_voov = abs(Hchk_voov(ub,ob,ob,ub) - Hex_voov(ub,ob,ob,ub));
    % 
    % fprintf('Error in H2A(amie) = %4.10f\n',sum(errA_voov(:)))
    % fprintf('Error in H2C(amie) = %4.10f\n',sum(errC_voov(:)))
    % fprintf('Error in H2B(amie) = %4.10f\n',sum(errB_voov(:)))
    % fprintf('Error in H2B(maie) = %4.10f\n',sum(errB_ovov(:)))
    % fprintf('Error in H2B(amei) = %4.10f\n',sum(errB_vovo(:)))
    % fprintf('Error in H2B(maei) = %4.10f\n',sum(errB_ovvo(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % 
    % % H(amef)
    % fprintf('\n==================== H(vovv) ====================\n')
    % 
    % % exact H(amef) in spinorbital form
    % Hex_vovv = HBar{2}{2,1,2,2};
    % Hex_ovvv = HBar{2}{1,2,2,2};
    % 
    % % spinorbital container
    % Hchk_vovv = zeros(size(Hex_vovv));
    % Hchk_ovvv = zeros(size(Hex_ovvv));
    % 
    % Hchk_vovv(ua,oa,ua,ua) = H2A.vovv;
    % Hchk_ovvv(oa,ua,ua,ua) = H2A.ovvv;
    % 
    % Hchk_vovv(ub,ob,ub,ub) = H2C.vovv;
    % Hchk_ovvv(ob,ub,ub,ub) = H2C.ovvv;
    % 
    % Hchk_vovv(ua,ob,ua,ub) = H2B.vovv;
    % Hchk_ovvv(oa,ub,ua,ub) = H2B.ovvv;
    % 
    % errA_vovv = abs(Hchk_vovv(ua,oa,ua,ua) - Hex_vovv(ua,oa,ua,ua));
    % errC_vovv = abs(Hchk_vovv(ub,ob,ub,ub) - Hex_vovv(ub,ob,ub,ub));
    % errB_vovv = abs(Hchk_vovv(ua,ob,ua,ub) - Hex_vovv(ua,ob,ua,ub));
    % 
    % errA_ovvv = abs(Hchk_ovvv(oa,ua,ua,ua) - Hex_ovvv(oa,ua,ua,ua));
    % errC_ovvv = abs(Hchk_ovvv(ob,ub,ub,ub) - Hex_ovvv(ob,ub,ub,ub));
    % errB_ovvv = abs(Hchk_ovvv(oa,ub,ua,ub) - Hex_ovvv(oa,ub,ua,ub));
    % 
    % fprintf('Error in H2A(amef) = %4.10f\n',sum(errA_vovv(:)))
    % fprintf('Error in H2C(amef) = %4.10f\n',sum(errC_vovv(:)))
    % fprintf('Error in H2B(amef) = %4.10f\n',sum(errB_vovv(:)))
    % fprintf('\n')
    % fprintf('Error in H2A(maef) = %4.10f\n',sum(errA_ovvv(:)))
    % fprintf('Error in H2C(maef) = %4.10f\n',sum(errC_ovvv(:)))
    % fprintf('Error in H2B(mafe) = %4.10f\n',sum(errB_ovvv(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % % H(mnie)
    % fprintf('\n==================== H(ooov) ====================\n')
    % 
    % % exact H(mnie) in spinorbital form
    % Hex_ooov = HBar{2}{1,1,1,2};
    % Hex_oovo = HBar{2}{1,1,2,1};
    % 
    % % spinorbital container
    % Hchk_ooov = zeros(size(Hex_ooov));
    % Hchk_oovo = zeros(size(Hex_oovo));
    % 
    % Hchk_ooov(oa,oa,oa,ua) = H2A.ooov;
    % Hchk_oovo(oa,oa,ua,oa) = H2A.oovo;
    % 
    % Hchk_ooov(ob,ob,ob,ub) = H2C.ooov;
    % Hchk_oovo(ob,ob,ub,ob) = H2C.oovo;
    % 
    % Hchk_ooov(oa,ob,oa,ub) = H2B.ooov;
    % Hchk_oovo(oa,ob,ua,ob) = H2B.oovo;
    % 
    % errA_ooov = abs(Hchk_ooov(oa,oa,oa,ua) - Hex_ooov(oa,oa,oa,ua));
    % errC_ooov = abs(Hchk_ooov(ob,ob,ob,ub) - Hex_ooov(ob,ob,ob,ub));
    % errB_ooov = abs(Hchk_ooov(oa,ob,oa,ub) - Hex_ooov(oa,ob,oa,ub));
    % 
    % errA_oovo = abs(Hchk_oovo(oa,oa,ua,oa) - Hex_oovo(oa,oa,ua,oa));
    % errC_oovo = abs(Hchk_oovo(ob,ob,ub,ob) - Hex_oovo(ob,ob,ub,ob));
    % errB_oovo = abs(Hchk_oovo(oa,ob,ua,ob) - Hex_oovo(oa,ob,ua,ob));
    % 
    % fprintf('Error in H2A(mnie) = %4.10f\n',sum(errA_ooov(:)))
    % fprintf('Error in H2C(mnie) = %4.10f\n',sum(errC_ooov(:)))
    % fprintf('Error in H2B(mnie) = %4.10f\n',sum(errB_ooov(:)))
    % fprintf('\n')
    % fprintf('Error in H2A(nmie) = %4.10f\n',sum(errA_oovo(:)))
    % fprintf('Error in H2C(nmie) = %4.10f\n',sum(errC_oovo(:)))
    % fprintf('Error in H2B(nmei) = %4.10f\n',sum(errB_oovo(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % fprintf('\n======================================== ++ 3-Body HBar ++ ========================================\n')
    % 
    % % H(abmije)
    % fprintf('\n=================== H(vvooov) ===================\n')
    % 
    % % exact H(abmije) in spinorbital form
    % Hex_vvooov = HBar{3}{2,2,1,1,1,2};
    % 
    % % spinorbital container
    % Hchk_vvooov = zeros(size(Hex_vvooov));
    % 
    % Hchk_vvooov(ua,ua,oa,oa,oa,ua) = H3A.vvooov;
    % Hchk_vvooov(ua,ua,ob,oa,oa,ub) = H3B.vvooov;
    % Hchk_vvooov(ua,ub,ob,oa,ob,ub) = H3C.vvooov;
    % Hchk_vvooov(ub,ub,ob,ob,ob,ub) = H3D.vvooov;
    % 
    % errA_vvooov = abs(Hchk_vvooov(ua,ua,oa,oa,oa,ua) - Hex_vvooov(ua,ua,oa,oa,oa,ua));
    % errB_vvooov = abs(Hchk_vvooov(ua,ua,ob,oa,oa,ub) - Hex_vvooov(ua,ua,ob,oa,oa,ub));
    % errC_vvooov = abs(Hchk_vvooov(ua,ub,ob,oa,ob,ub) - Hex_vvooov(ua,ub,ob,oa,ob,ub));
    % errD_vvooov = abs(Hchk_vvooov(ub,ub,ob,ob,ob,ub) - Hex_vvooov(ub,ub,ob,ob,ob,ub));
    % 
    % fprintf('Error in H3A(abmije) = %4.10f\n',sum(errA_vvooov(:)))
    % fprintf('Error in H3B(abmije) = %4.10f\n',sum(errB_vvooov(:)))
    % fprintf('Error in H3C(abmije) = %4.10f\n',sum(errC_vvooov(:)))
    % fprintf('Error in H3D(abmije) = %4.10f\n',sum(errD_vvooov(:)))
    % 
    % fprintf('\n')
    % 
    % % exact H(ambiej) in spinorbital form
    % Hex_vovovo = HBar{3}{2,1,2,1,2,1};
    % 
    % % spinorbital container
    % Hchk_vovovo = zeros(size(Hex_vovovo));
    % 
    % Hchk_vovovo(ua,oa,ua,oa,ua,oa) = H3A.vovovo;
    % Hchk_vovovo(ua,oa,ub,oa,ua,ob) = H3B.vovovo;
    % Hchk_vovovo(ua,ob,ub,oa,ub,ob) = H3C.vovovo;
    % Hchk_vovovo(ub,ob,ub,ob,ub,ob) = H3D.vovovo;
    % 
    % errA_vovovo = abs(Hchk_vovovo(ua,oa,ua,oa,ua,oa) - Hex_vovovo(ua,oa,ua,oa,ua,oa));
    % errB_vovovo = abs(Hchk_vovovo(ua,oa,ub,oa,ua,ob) - Hex_vovovo(ua,oa,ub,oa,ua,ob));
    % errC_vovovo = abs(Hchk_vovovo(ua,ob,ub,oa,ub,ob) - Hex_vovovo(ua,ob,ub,oa,ub,ob));
    % errD_vovovo = abs(Hchk_vovovo(ub,ob,ub,ob,ub,ob) - Hex_vovovo(ub,ob,ub,ob,ub,ob));
    % 
    % fprintf('Error in H3A(ambiej) = %4.10f\n',sum(errA_vovovo(:)))
    % fprintf('Error in H3B(ambiej) = %4.10f\n',sum(errB_vovovo(:)))
    % fprintf('Error in H3C(ambiej) = %4.10f\n',sum(errC_vovovo(:)))
    % fprintf('Error in H3D(ambiej) = %4.10f\n',sum(errD_vovovo(:)))
    % 
    % fprintf('\n')
    % 
    % % exact H(mabeij) in spinorbital form
    % Hex_ovvvoo = HBar{3}{1,2,2,2,1,1};
    % 
    % % spinorbital container
    % Hchk_ovvvoo = zeros(size(Hex_ovvvoo));
    % 
    % Hchk_ovvvoo(oa,ua,ua,ua,oa,oa) = H3A.ovvvoo;
    % Hchk_ovvvoo(oa,ua,ub,ua,oa,ob) = H3B.ovvvoo;
    % Hchk_ovvvoo(oa,ub,ub,ua,ob,ob) = H3C.ovvvoo;
    % Hchk_ovvvoo(ob,ub,ub,ub,ob,ob) = H3D.ovvvoo;
    % 
    % errA_ovvvoo = abs(Hchk_ovvvoo(oa,ua,ua,ua,oa,oa) - Hex_ovvvoo(oa,ua,ua,ua,oa,oa));
    % errB_ovvvoo = abs(Hchk_ovvvoo(oa,ua,ub,ua,oa,ob) - Hex_ovvvoo(oa,ua,ub,ua,oa,ob));
    % errC_ovvvoo = abs(Hchk_ovvvoo(oa,ub,ub,ua,ob,ob) - Hex_ovvvoo(oa,ub,ub,ua,ob,ob));
    % errD_ovvvoo = abs(Hchk_ovvvoo(ob,ub,ub,ub,ob,ob) - Hex_ovvvoo(ob,ub,ub,ub,ob,ob));
    % 
    % fprintf('Error in H3A(mabeij) = %4.10f\n',sum(errA_ovvvoo(:)))
    % fprintf('Error in H3B(mabeij) = %4.10f\n',sum(errB_ovvvoo(:)))
    % fprintf('Error in H3C(mabeij) = %4.10f\n',sum(errC_ovvvoo(:)))
    % fprintf('Error in H3D(mabeij) = %4.10f\n',sum(errD_ovvvoo(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % % H(abmief)
    % fprintf('\n=================== H(vvoovv) ===================\n')
    % 
    % % exact H(abmief) in spinorbital form
    % Hex_vvoovv = HBar{3}{2,2,1,1,2,2};
    % 
    % % spinorbital container
    % Hchk_vvoovv = zeros(size(Hex_vvoovv));
    % 
    % Hchk_vvoovv(ua,ua,oa,oa,ua,ua) = H3A.vvoovv;
    % Hchk_vvoovv(ua,ua,ob,oa,ua,ub) = H3B.vvoovv;
    % Hchk_vvoovv(ua,ub,ob,oa,ub,ub) = H3C.vvoovv;
    % Hchk_vvoovv(ub,ub,ob,ob,ub,ub) = H3D.vvoovv;
    % 
    % errA_vvoovv = abs(Hchk_vvoovv(ua,ua,oa,oa,ua,ua) - Hex_vvoovv(ua,ua,oa,oa,ua,ua));
    % errB_vvoovv = abs(Hchk_vvoovv(ua,ua,ob,oa,ua,ub) - Hex_vvoovv(ua,ua,ob,oa,ua,ub));
    % errC_vvoovv = abs(Hchk_vvoovv(ua,ub,ob,oa,ub,ub) - Hex_vvoovv(ua,ub,ob,oa,ub,ub));
    % errD_vvoovv = abs(Hchk_vvoovv(ub,ub,ob,ob,ub,ub) - Hex_vvoovv(ub,ub,ob,ob,ub,ub));
    % 
    % fprintf('Error in H3A(abmief) = %4.10f\n',sum(errA_vvoovv(:)))
    % fprintf('Error in H3B(abmief) = %4.10f\n',sum(errB_vvoovv(:)))
    % fprintf('Error in H3C(abmief) = %4.10f\n',sum(errC_vvoovv(:)))
    % fprintf('Error in H3D(abmief) = %4.10f\n',sum(errD_vvoovv(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % % H(amnije)
    % fprintf('\n=================== H(voooov) ===================\n')
    % 
    % % exact H(amnije) in spinorbital form
    % Hex_voooov = HBar{3}{2,1,1,1,1,2};
    % 
    % % spinorbital container
    % Hchk_voooov = zeros(size(Hex_voooov));
    % 
    % Hchk_voooov(ua,oa,oa,oa,oa,ua) = H3A.voooov;
    % Hchk_voooov(ua,oa,ob,oa,oa,ub) = H3B.voooov;
    % Hchk_voooov(ua,ob,ob,oa,ob,ub) = H3C.voooov;
    % Hchk_voooov(ub,ob,ob,ob,ob,ub) = H3D.voooov;
    % 
    % errA_voooov = abs(Hchk_voooov(ua,oa,oa,oa,oa,ua) - Hex_voooov(ua,oa,oa,oa,oa,ua));
    % errB_voooov = abs(Hchk_voooov(ua,oa,ob,oa,oa,ub) - Hex_voooov(ua,oa,ob,oa,oa,ub));
    % errC_voooov = abs(Hchk_voooov(ua,ob,ob,oa,ob,ub) - Hex_voooov(ua,ob,ob,oa,ob,ub));
    % errD_voooov = abs(Hchk_voooov(ub,ob,ob,ob,ob,ub) - Hex_voooov(ub,ob,ob,ob,ob,ub));
    % 
    % fprintf('Error in H3A(amnije) = %4.10f\n',sum(errA_voooov(:)))
    % fprintf('Error in H3B(amnije) = %4.10f\n',sum(errB_voooov(:)))
    % fprintf('Error in H3C(amnije) = %4.10f\n',sum(errC_voooov(:)))
    % fprintf('Error in H3D(amnije) = %4.10f\n',sum(errD_voooov(:)))
    % 
    % fprintf('=================================================\n')
    % 
    % 
    % 
    % %%


