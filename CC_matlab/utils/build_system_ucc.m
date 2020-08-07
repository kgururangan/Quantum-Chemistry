function [sys] = build_system_ucc(e1int,e2int,Vnuc,Nocc_alpha,Nocc_beta,nfzc)

    % THIS IS ONLY SET UP FOR RHF AT THE MOMENT. e2int AND e1int HAVE
    % DIMENSIONS OF Norb AND Norb. FOR A UHF REFERENCE, WE WOULD HAVE 
    % AN e1int AND e2int FOR ALPHA AND BETA ORBITALS, SEPARATELY.
    
    % i.e. the e1int would be of size 2*Norb x 2*Norb and e2int would be of
    % size 2*Norb x 2*Norb x 2*Norb x 2*Norb and we would extract
    
    % e1intA = e1int(1:2:2*Norb,1:2:2*Norb)
    % e1intB = e1int(2:2:2*Norb,2:2:2*Norb)
    % e2intA = e2int(1:2:2*Norb,1:2:2*Norb,1:2:2*Norb,1:2:2*Norb)
    % e2intB = e2int(2:2:2*Norb,2:2:2*Norb,2:2:2*Norb,2:2:2*Norb)
    % vA = e2intA - permute(e2intA,[1,2,4,3])
    % vB = e2int(1:2:2*Norb,2:2:2*Norb,1:2:2*Norb,2:2:2*Norb)
    % vC = e2intB - permute(e2intB,[1,2,4,3])
        
    % total number of spatial orbitals in system
    Norb = size(e1int,1);

    % Calculate the scf energy
    hf_energy = 0.0;
    for i = 1:Nocc_alpha
        hf_energy = hf_energy + e1int(i,i);
    end
    for i = 1:Nocc_beta
        hf_energy = hf_energy + e1int(i,i);
    end
    for i = 1:Nocc_alpha
        for j = 1:Nocc_beta
            hf_energy = hf_energy + e2int(i,j,i,j);
        end
    end
    for i = 1:Nocc_alpha
        for j = 1:Nocc_alpha
            hf_energy = hf_energy + 0.5*(e2int(i,j,i,j)-e2int(i,j,j,i));
        end
    end
    for i = 1:Nocc_beta
        for j = 1:Nocc_beta
            hf_energy = hf_energy + 0.5*(e2int(i,j,i,j)-e2int(i,j,j,i));
        end
    end
    sys.Escf = hf_energy + Vnuc;
    
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

    % store full unfrozen spatial orbital integrals (for Goldstone stuff)
    sys.e1int = e1int;
    sys.e2int = e2int;
    sys.Fock_a = FA;
    sys.Fock_b = FB;
    sys.vA = VA;
    sys.vB = VB;
    sys.vC = VC;
    
    % freeze orbitals and define slicing vectors
    iocc_alpha = nfzc+1:Nocc_alpha;
    iocc_beta = nfzc+1:Nocc_beta;
    ivir_alpha = Nocc_alpha+1:Norb;
    ivir_beta = Nocc_beta+1:Norb;

    Nocc_alpha = length(iocc_alpha); Nocc_beta = length(iocc_beta);
    Nvir_alpha = length(ivir_alpha); Nvir_beta = length(ivir_beta);
    
    % number of correlated electrons    
    Nelec = Nocc_alpha + Nocc_beta;

    % populate system structure for spin-integrated post-HF methods
    sys.Norb = Norb;
    sys.Nelec = Nelec;
    sys.Nocc_alpha = length(iocc_alpha);
    sys.Nocc_beta = length(iocc_beta);
    sys.Nvir_alpha = length(ivir_alpha);
    sys.Nvir_beta = length(ivir_beta);
    sys.iocc_alpha = iocc_alpha;
    sys.ivir_alpha = ivir_alpha;
    sys.iocc_beta = iocc_beta;
    sys.ivir_beta = ivir_beta;

    sys.singles_dim = Nocc_alpha*Nvir_alpha + ...
                      Nocc_beta*Nvir_beta;

    sys.doubles_dim = Nocc_alpha*Nvir_alpha + ...
                      Nocc_beta*Nvir_beta + ...
                      Nocc_alpha^2*Nvir_alpha^2 + ...
                      Nocc_alpha*Nocc_beta*Nvir_alpha*Nvir_beta + ...
                      Nocc_beta^2*Nvir_beta^2;

    sys.triples_dim = Nocc_alpha*Nvir_alpha + ...
                      Nocc_beta*Nvir_beta + ...
                      Nocc_alpha^2*Nvir_alpha^2 + ...
                      Nocc_alpha*Nocc_beta*Nvir_alpha*Nvir_beta + ...
                      Nocc_beta^2*Nvir_beta^2 + ...
                      Nocc_alpha^3*Nvir_alpha^3 + ...
                      Nocc_alpha^2*Nocc_beta*Nvir_alpha^2*Nvir_beta + ...
                      Nocc_alpha*Nocc_beta^2*Nvir_alpha*Nvir_beta^2 + ...
                      Nocc_beta^3*Nvir_beta^3;

    post1a = 1:Nocc_alpha*Nvir_alpha;
    post1b = [post1a(end)+1:post1a(end)+Nocc_beta*Nvir_beta];
    post2a = [post1b(end)+1:post1b(end)+Nocc_alpha^2*Nvir_beta^2];
    post2b = [post2a(end)+1:post2a(end)+Nocc_alpha*Nocc_beta*Nvir_alpha*Nvir_beta];
    post2c = [post2b(end)+1:post2b(end)+Nocc_beta^2*Nvir_beta^2];
    sys.posv = {post1a, post1b, post2a, post2b, post2c};

    sys.size = {[Nvir_alpha, Nocc_alpha], [Nvir_beta, Nocc_beta], ...
                [Nvir_alpha, Nvir_alpha, Nocc_alpha, Nocc_alpha],...
                [Nvir_alpha, Nvir_beta, Nocc_alpha, Nocc_beta],...
                [Nvir_beta, Nvir_beta, Nocc_beta, Nocc_beta]};

    % diagonal fock masks
    Zocc_alpha = ones(sys.Nocc_alpha) - eye(sys.Nocc_alpha);
    Zvir_alpha = ones(sys.Nvir_alpha) - eye(sys.Nvir_alpha);
    Zocc_beta = ones(sys.Nocc_beta) - eye(sys.Nocc_beta);
    Zvir_beta = ones(sys.Nvir_beta) - eye(sys.Nvir_beta);
    
    % f1A
    sys.fa_oo = FA(iocc_alpha, iocc_alpha);  sys.fa_oo_masked = sys.fa_oo.*Zocc_alpha;
    sys.fa_vv = FA(ivir_alpha, ivir_alpha);  sys.fa_vv_masked = sys.fa_vv.*Zvir_alpha;
    sys.fa_ov = FA(iocc_alpha, ivir_alpha);
    sys.fa_vo = FA(ivir_alpha, iocc_alpha);

    % f1B
    sys.fb_oo = FB(iocc_beta, iocc_beta);   sys.fb_oo_masked = sys.fb_oo.*Zocc_beta;
    sys.fb_vv = FB(ivir_beta, ivir_beta);   sys.fb_vv_masked = sys.fb_vv.*Zvir_beta;
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

    % Goldstone e2int slicing
    sys.e2int_voov = e2int(ivir_alpha,iocc_alpha,iocc_alpha,ivir_alpha);
    sys.e2int_vovo = e2int(ivir_alpha,iocc_alpha,ivir_alpha,iocc_alpha);

end