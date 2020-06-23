function [sys] = build_system_ucc(e1int,e2int,Nocc_alpha,Nocc_beta)

        
    Norb = size(e1int,1);
    
    Nelec = Nocc_alpha + Nocc_beta;
    
    iocc_alpha = 1:Nocc_alpha;
    iocc_beta = 1:Nocc_beta;
    ivir_alpha = Nocc_alpha+1:Norb;
    ivir_beta = Nocc_beta+1:Norb;

    % populate system structure for spin-integrated post-HF methods
    sys.Nelec = Nelec;
    sys.Nocc_alpha = length(iocc_alpha);
    sys.Nocc_beta = length(iocc_beta);
    sys.Nvir_alpha = length(ivir_alpha);
    sys.Nvir_beta = length(ivir_beta);
    sys.iocc_alpha = iocc_alpha;
    sys.ivir_alpha = ivir_alpha;
    sys.iocc_beta = iocc_beta;
    sys.ivir_beta = ivir_beta;
    
    % build antisymmetrized twobody integrals
    VA = e2int - permute(e2int,[1,2,4,3]);
    VB = e2int;
    VC = e2int - permute(e2int,[1,2,4,3]);
    
    % build Fock matrices using F = Z + G for different spin cases
    FA = zeros(Norb); FB = zeros(Norb);
    for p = 1:Norb
        for q = 1:Norb
               FA(p,q) = e1int(p,q);
               FB(p,q) = e1int(p,q);
               for i = 1:Nocc_beta
                   FA(p,q) = FA(p,q) + VA(p,i,q,i) + VB(p,i,q,i);
               end
               for i = 1:Nocc_alpha
                   FB(p,q) = FB(p,q) + VC(p,i,q,i) + VB(i,p,i,q);
               end
        end
    end
    
    
    % f1A
    sys.fa_oo = FA(iocc_alpha, iocc_alpha); 
    sys.fa_vv = FA(ivir_alpha, ivir_alpha);
    sys.fa_ov = FA(iocc_alpha, ivir_alpha);
    sys.fa_vo = FA(ivir_alpha, iocc_alpha);

    % f1B
    sys.fb_oo = FB(iocc_beta, iocc_beta);
    sys.fb_vv = FB(ivir_beta, ivir_beta);
    sys.fb_ov = FB(iocc_beta, ivir_beta);
    sys.fb_vo = FB(ivir_beta, iocc_beta);

    % v2A
    sys.vA_oooo = VA(iocc_alpha, iocc_alpha, iocc_alpha, iocc_alpha);
    sys.vA_ooov = VA(iocc_alpha, iocc_alpha, iocc_alpha, ivir_alpha);
    sys.vA_oovo = VA(iocc_alpha, iocc_alpha, ivir_alpha, iocc_alpha);
    sys.vA_ovoo = VA(iocc_alpha, ivir_alpha, iocc_alpha, iocc_alpha);
    sys.vA_vooo = VA(ivir_alpha, iocc_alpha, iocc_alpha, iocc_alpha);
    sys.vA_oovv = VA(iocc_alpha, iocc_alpha, ivir_alpha, ivir_alpha);
    sys.vA_ovov = VA(iocc_alpha, ivir_alpha, iocc_alpha, ivir_alpha);
    sys.vA_voov = VA(ivir_alpha, iocc_alpha, iocc_alpha, ivir_alpha);
    sys.vA_vovo = VA(ivir_alpha, iocc_alpha, ivir_alpha, iocc_alpha);
    sys.vA_vvoo = VA(ivir_alpha, ivir_alpha, iocc_alpha, iocc_alpha);
    sys.vA_ovvo = VA(iocc_alpha, ivir_alpha, ivir_alpha, iocc_alpha);
    sys.vA_ovvv = VA(iocc_alpha, ivir_alpha, ivir_alpha, ivir_alpha);
    sys.vA_vovv = VA(ivir_alpha, iocc_alpha, ivir_alpha, ivir_alpha);
    sys.vA_vvov = VA(ivir_alpha, ivir_alpha, iocc_alpha, ivir_alpha);
    sys.vA_vvvo = VA(ivir_alpha, ivir_alpha, ivir_alpha, iocc_alpha);
    sys.vA_vvvv = VA(ivir_alpha, ivir_alpha, ivir_alpha, ivir_alpha);

    % v2B
    sys.vB_oooo = VB(iocc_alpha, iocc_beta, iocc_alpha, iocc_beta);
    sys.vB_ooov = VB(iocc_alpha, iocc_beta, iocc_alpha, ivir_beta);
    sys.vB_oovo = VB(iocc_alpha, iocc_beta, ivir_alpha, iocc_beta);
    sys.vB_ovoo = VB(iocc_alpha, ivir_beta, iocc_alpha, iocc_beta);
    sys.vB_vooo = VB(ivir_alpha, iocc_beta, iocc_alpha, iocc_beta);
    sys.vB_oovv = VB(iocc_alpha, iocc_beta, ivir_alpha, ivir_beta);
    sys.vB_ovov = VB(iocc_alpha, ivir_beta, iocc_alpha, ivir_beta);
    sys.vB_voov = VB(ivir_alpha, iocc_beta, iocc_alpha, ivir_beta);
    sys.vB_vovo = VB(ivir_alpha, iocc_beta, ivir_alpha, iocc_beta);
    sys.vB_vvoo = VB(ivir_alpha, ivir_beta, iocc_alpha, iocc_beta);
    sys.vB_ovvo = VB(iocc_alpha, ivir_beta, ivir_alpha, iocc_beta);
    sys.vB_ovvv = VB(iocc_alpha, ivir_beta, ivir_alpha, ivir_beta);
    sys.vB_vovv = VB(ivir_alpha, iocc_beta, ivir_alpha, ivir_beta);
    sys.vB_vvov = VB(ivir_alpha, ivir_beta, iocc_alpha, ivir_beta);
    sys.vB_vvvo = VB(ivir_alpha, ivir_beta, ivir_alpha, iocc_beta);
    sys.vB_vvvv = VB(ivir_alpha, ivir_beta, ivir_alpha, ivir_beta);

    % v2C
    sys.vC_oooo = VC(iocc_beta, iocc_beta, iocc_beta, iocc_beta);
    sys.vC_ooov = VC(iocc_beta, iocc_beta, iocc_beta, ivir_beta);
    sys.vC_oovo = VC(iocc_beta, iocc_beta, ivir_beta, iocc_beta);
    sys.vC_ovoo = VC(iocc_beta, ivir_beta, iocc_beta, iocc_beta);
    sys.vC_vooo = VC(ivir_beta, iocc_beta, iocc_beta, iocc_beta);
    sys.vC_oovv = VC(iocc_beta, iocc_beta, ivir_beta, ivir_beta);
    sys.vC_ovov = VC(iocc_beta, ivir_beta, iocc_beta, ivir_beta);
    sys.vC_voov = VC(ivir_beta, iocc_beta, iocc_beta, ivir_beta);
    sys.vC_vovo = VC(ivir_beta, iocc_beta, ivir_beta, iocc_beta);
    sys.vC_vvoo = VC(ivir_beta, ivir_beta, iocc_beta, iocc_beta);
    sys.vC_ovvo = VC(iocc_beta, ivir_beta, ivir_beta, iocc_beta);
    sys.vC_ovvv = VC(iocc_beta, ivir_beta, ivir_beta, ivir_beta);
    sys.vC_vovv = VC(ivir_beta, iocc_beta, ivir_beta, ivir_beta);
    sys.vC_vvov = VC(ivir_beta, ivir_beta, iocc_beta, ivir_beta);
    sys.vC_vvvo = VC(ivir_beta, ivir_beta, ivir_beta, iocc_beta);
    sys.vC_vvvv = VC(ivir_beta, ivir_beta, ivir_beta, ivir_beta);

end