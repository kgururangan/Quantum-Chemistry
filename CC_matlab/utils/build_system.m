function [sys] = build_system(VM,FM,occ,unocc,Nelec)

% This function is used to convert spinorbital electron-repulsion and Fock
% matrices into the system struct used in UCC codes

    if occ(1) == 0
        occ = occ + 1;
        unocc = unocc + 1;
    end
    

    ialpha = [occ(1:2:end-1),unocc(1:2:end-1)];
    ibeta = [occ(2:2:end),unocc(2:2:end)];

    iocc_alpha = occ(1:2:end-1);
    iocc_beta = occ(2:2:end);
    ivir_alpha = unocc(1:2:end);
    ivir_beta = unocc(2:2:end);

    Nocc_alpha = length(iocc_alpha); Nocc_beta = length(iocc_beta);
    Nvir_alpha = length(ivir_alpha); Nvir_beta = length(ivir_beta);

    % populate system structure for spin-integrated post-HF methods
    sys.Nelec = Nelec;
    sys.Nocc = length(occ)/2;
    sys.Nvir = length(unocc)/2;
    sys.Nov = sys.Nocc*sys.Nvir;
    sys.Nocc_alpha = Nocc_alpha;
    sys.Nvir_alpha = Nvir_alpha;
    sys.Nocc_beta = Nocc_beta;
    sys.Nvir_beta = Nvir_beta;
    sys.ialpha = ialpha;
    sys.ibeta = ibeta;
    sys.iocc = occ;
    sys.ivir = unocc;
    sys.iocc_alpha = iocc_alpha;
    sys.iocc_beta = iocc_beta;
    sys.ivir_alpha = ivir_alpha;
    sys.ivir_beta = ivir_beta;
    sys.VM = VM;
    sys.FM = FM;

    % f1A
    sys.fa_oo = FM(iocc_alpha, iocc_alpha); 
    sys.fa_vv = FM(ivir_alpha, ivir_alpha);
    sys.fa_ov = FM(iocc_alpha, ivir_alpha);
    sys.fa_vo = FM(ivir_alpha, iocc_alpha);

    % f1B
    sys.fb_oo = FM(iocc_beta, iocc_beta);
    sys.fb_vv = FM(ivir_beta, ivir_beta);
    sys.fb_ov = FM(iocc_beta, ivir_beta);
    sys.fb_vo = FM(ivir_beta, iocc_beta);

    % v2A
    sys.vA_oooo = VM(iocc_alpha, iocc_alpha, iocc_alpha, iocc_alpha);
    sys.vA_ooov = VM(iocc_alpha, iocc_alpha, iocc_alpha, ivir_alpha);
    sys.vA_oovo = VM(iocc_alpha, iocc_alpha, ivir_alpha, iocc_alpha);
    sys.vA_ovoo = VM(iocc_alpha, ivir_alpha, iocc_alpha, iocc_alpha);
    sys.vA_vooo = VM(ivir_alpha, iocc_alpha, iocc_alpha, iocc_alpha);
    sys.vA_oovv = VM(iocc_alpha, iocc_alpha, ivir_alpha, ivir_alpha);
    sys.vA_ovov = VM(iocc_alpha, ivir_alpha, iocc_alpha, ivir_alpha);
    sys.vA_voov = VM(ivir_alpha, iocc_alpha, iocc_alpha, ivir_alpha);
    sys.vA_vovo = VM(ivir_alpha, iocc_alpha, ivir_alpha, iocc_alpha);
    sys.vA_vvoo = VM(ivir_alpha, ivir_alpha, iocc_alpha, iocc_alpha);
    sys.vA_ovvo = VM(iocc_alpha, ivir_alpha, ivir_alpha, iocc_alpha);
    sys.vA_ovvv = VM(iocc_alpha, ivir_alpha, ivir_alpha, ivir_alpha);
    sys.vA_vovv = VM(ivir_alpha, iocc_alpha, ivir_alpha, ivir_alpha);
    sys.vA_vvov = VM(ivir_alpha, ivir_alpha, iocc_alpha, ivir_alpha);
    sys.vA_vvvo = VM(ivir_alpha, ivir_alpha, ivir_alpha, iocc_alpha);
    sys.vA_vvvv = VM(ivir_alpha, ivir_alpha, ivir_alpha, ivir_alpha);

    % v2B
    sys.vB_oooo = VM(iocc_alpha, iocc_beta, iocc_alpha, iocc_beta);
    sys.vB_ooov = VM(iocc_alpha, iocc_beta, iocc_alpha, ivir_beta);
    sys.vB_oovo = VM(iocc_alpha, iocc_beta, ivir_alpha, iocc_beta);
    sys.vB_ovoo = VM(iocc_alpha, ivir_beta, iocc_alpha, iocc_beta);
    sys.vB_vooo = VM(ivir_alpha, iocc_beta, iocc_alpha, iocc_beta);
    sys.vB_oovv = VM(iocc_alpha, iocc_beta, ivir_alpha, ivir_beta);
    sys.vB_ovov = VM(iocc_alpha, ivir_beta, iocc_alpha, ivir_beta);
    sys.vB_voov = VM(ivir_alpha, iocc_beta, iocc_alpha, ivir_beta);
    sys.vB_vovo = VM(ivir_alpha, iocc_beta, ivir_alpha, iocc_beta);
    sys.vB_vvoo = VM(ivir_alpha, ivir_beta, iocc_alpha, iocc_beta);
    sys.vB_ovvo = VM(iocc_alpha, ivir_beta, ivir_alpha, iocc_beta);
    sys.vB_ovvv = VM(iocc_alpha, ivir_beta, ivir_alpha, ivir_beta);
    sys.vB_vovv = VM(ivir_alpha, iocc_beta, ivir_alpha, ivir_beta);
    sys.vB_vvov = VM(ivir_alpha, ivir_beta, iocc_alpha, ivir_beta);
    sys.vB_vvvo = VM(ivir_alpha, ivir_beta, ivir_alpha, iocc_beta);
    sys.vB_vvvv = VM(ivir_alpha, ivir_beta, ivir_alpha, ivir_beta);

    % v2C
    sys.vC_oooo = VM(iocc_beta, iocc_beta, iocc_beta, iocc_beta);
    sys.vC_ooov = VM(iocc_beta, iocc_beta, iocc_beta, ivir_beta);
    sys.vC_oovo = VM(iocc_beta, iocc_beta, ivir_beta, iocc_beta);
    sys.vC_ovoo = VM(iocc_beta, ivir_beta, iocc_beta, iocc_beta);
    sys.vC_vooo = VM(ivir_beta, iocc_beta, iocc_beta, iocc_beta);
    sys.vC_oovv = VM(iocc_beta, iocc_beta, ivir_beta, ivir_beta);
    sys.vC_ovov = VM(iocc_beta, ivir_beta, iocc_beta, ivir_beta);
    sys.vC_voov = VM(ivir_beta, iocc_beta, iocc_beta, ivir_beta);
    sys.vC_vovo = VM(ivir_beta, iocc_beta, ivir_beta, iocc_beta);
    sys.vC_vvoo = VM(ivir_beta, ivir_beta, iocc_beta, iocc_beta);
    sys.vC_ovvo = VM(iocc_beta, ivir_beta, ivir_beta, iocc_beta);
    sys.vC_ovvv = VM(iocc_beta, ivir_beta, ivir_beta, ivir_beta);
    sys.vC_vovv = VM(ivir_beta, iocc_beta, ivir_beta, ivir_beta);
    sys.vC_vvov = VM(ivir_beta, ivir_beta, iocc_beta, ivir_beta);
    sys.vC_vvvo = VM(ivir_beta, ivir_beta, ivir_beta, iocc_beta);
    sys.vC_vvvv = VM(ivir_beta, ivir_beta, ivir_beta, ivir_beta);

end

