function [sys] = build_system_ucc(e1int,e2int,Vnuc,Nocc_alpha,Nocc_beta,nfzc,nact_h_alpha,nact_p_alpha,nact_h_beta,nact_p_beta,flag_act_scheme)

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

    if ~exist('nact_h_alpha','var') 
        nact_h_alpha = 10000;
    end
    if ~exist('nact_h_beta','var') 
        nact_h_beta = 10000;
    end
    if ~exist('nact_p_alpha','var')
        nact_p_alpha = 10000;
    end
    if ~exist('nact_p_beta','var')
        nact_p_beta = 10000;
    end
    if ~exist('flag_act_scheme','var')
        flag_act_scheme = 0;
    end
        
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
    
    % freeze orbitals and define V, F slicing vectors
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
    
    % active space
    Nact_h_alpha = min(Nocc_alpha,nact_h_alpha);
    Nact_p_alpha = min(Nvir_alpha,nact_p_alpha);
    Nact_h_beta = min(Nocc_beta,nact_h_beta);
    Nact_p_beta = min(Nvir_beta,nact_p_beta);
    Nunact_h_alpha = Nocc_alpha - Nact_h_alpha;
    Nunact_p_alpha = Nvir_alpha - Nact_p_alpha;
    Nunact_h_beta = Nocc_beta - Nact_h_beta;
    Nunact_p_beta = Nvir_beta - Nact_p_beta;
    
    % store active space dimensions
    sys.Nact_h_alpha = Nact_h_alpha;
    sys.Nact_h_beta = Nact_h_beta;
    sys.Nact_p_alpha = Nact_p_alpha;
    sys.Nact_p_beta = Nact_p_beta;
    sys.Nunact_h_alpha = Nunact_h_alpha;
    sys.Nunact_h_beta = Nunact_h_beta;
    sys.Nunact_p_alpha = Nunact_p_alpha;
    sys.Nunact_p_beta = Nunact_p_beta;
   
    % VA, VB, VC active space slicing vectors
    iact_p_alpha = [Nocc_alpha+nfzc+1 : Nocc_alpha+nfzc+Nact_p_alpha];
    iunact_p_alpha = [iact_p_alpha(end)+1:Norb];
    iact_h_alpha = [Nocc_alpha+nfzc-Nact_h_alpha+1 : Nocc_alpha+nfzc];
    iunact_h_alpha = [nfzc+1:iact_h_alpha(1)-1];
    iact_p_beta = [Nocc_beta+nfzc+1 : Nocc_beta+nfzc+Nact_p_beta];
    iunact_p_beta = [iact_p_beta(end)+1:Norb];
    iact_h_beta = [Nocc_beta+nfzc-Nact_h_beta+1 : Nocc_beta+nfzc];
    iunact_h_beta = [nfzc+1:iact_h_beta(1)-1];
   
    sys.iPA = iact_p_alpha;     sys.iPB = iact_p_beta;
    sys.ipA = iunact_p_alpha;   sys.ipB = iunact_p_beta;
    sys.iHA = iact_h_alpha;     sys.iHB = iact_h_beta;
    sys.ihA = iunact_h_alpha;   sys.ihB = iunact_h_beta;
    
    % T vector active slicing indices
    sys.HA = (Nocc_alpha-Nact_h_alpha+1):Nocc_alpha;
    sys.hA = 1:(iact_h_alpha(1)-nfzc-1);
    sys.PA = 1:Nact_p_alpha;
    sys.pA = (Nact_p_alpha+1):Nvir_alpha;
    sys.HB = (Nocc_beta-Nact_h_beta+1):Nocc_beta;
    sys.hB = 1:(iact_h_beta(1)-nfzc-1);
    sys.PB = 1:Nact_p_beta;
    sys.pB = (Nact_p_beta+1):Nvir_beta;
    
    % calculate Hilbert subspace dimension sizes
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
    
    switch flag_act_scheme
        
        case 0
            
            sys.size = {[Nvir_alpha, Nocc_alpha], [Nvir_beta, Nocc_beta], ...
                        [Nvir_alpha, Nvir_alpha, Nocc_alpha, Nocc_alpha],...
                        [Nvir_alpha, Nvir_beta, Nocc_alpha, Nocc_beta],...
                        [Nvir_beta, Nvir_beta, Nocc_beta, Nocc_beta],...
                        [Nvir_alpha, Nvir_alpha, Nvir_alpha, Nocc_alpha, Nocc_alpha, Nocc_alpha],...
                        [Nvir_alpha, Nvir_alpha, Nvir_beta, Nocc_alpha, Nocc_alpha, Nocc_beta],...
                        [Nvir_alpha, Nvir_beta, Nvir_beta, Nocc_alpha, Nocc_beta, Nocc_beta],...
                        [Nvir_beta, Nvir_beta, Nvir_beta, Nocc_beta, Nocc_beta, Nocc_beta]};
            
            if Norb < 40 % allocate full triples

                post3a = [post2c(end)+1:post2c(end)+Nocc_alpha^3*Nvir_alpha^3];

                post3b = [post3a(end)+1:post3a(end)+Nocc_alpha^2*Nocc_beta*Nvir_alpha^2*Nvir_beta];

                post3c = [post3b(end)+1:post3b(end)+Nocc_alpha*Nocc_beta^2*Nvir_alpha*Nvir_beta^2];

                post3d = [post3c(end)+1:post3c(end)+Nocc_beta^3*Nvir_beta^3];

                sys.posv = {post1a, post1b, post2a, post2b, post2c, post3a, post3b, post3c, post3d};

            else
                fprintf('Due to problem size, full triples space not allocated\n')
                sys.posv = {post1a, post1b, post2a, post2b, post2c};
            end
            
        case 2
            
                % IJK^ABC = IJK^ABC + 
                %           IJK~^ABC~ + 
                %           IJ~K~^AB~C~ + 
                %           I~J~K~^A~B~C~
                num_proj1 =  Nact_h_alpha^3*Nact_p_alpha^3 + ...
                             Nact_h_alpha^2*Nact_h_beta*Nact_p_alpha^2*Nact_p_beta + ...
                             Nact_h_alpha*Nact_h_beta^2*Nact_p_alpha*Nact_p_beta^2 + ...
                             Nact_h_beta^3*Nact_p_beta^3;
                
                % iJK^ABC = iJK^ABC +
                %           iJK~^ABC~ + IJk~^ABC~ + 
                %           iJ~K~^AB~C~ + Ij~K~^AB~C~ +
                %           i~J~K~^A~B~C~
                num_proj2 = Nunact_h_alpha*Nact_h_alpha^2*Nact_p_alpha^3 + ... % A
                          Nunact_h_alpha*Nact_h_alpha*Nact_h_beta*Nact_p_alpha^2*Nact_p_beta + ...
                          Nunact_h_beta*Nact_h_alpha^2*Nact_p_alpha^2*Nact_p_beta + ... % B
                          Nunact_h_alpha*Nact_h_beta^2*Nact_p_alpha*Nact_p_beta^2 + ...
                          Nact_h_alpha*Nact_h_beta*Nunact_h_beta*Nact_p_alpha*Nact_p_beta^2 + ... % C
                          Nunact_h_beta*Nact_h_beta^2*Nact_p_beta^3; % D

                % IJK^ABc = IJK^ABc +
                %           IJK~^ABc~ + IKJ~^AcB~ + 
                %           IJ~K~^AB~c~ + KI~J~^cA~B~ +
                %           I~J~K~^A~B~c~
                num_proj3 = Nact_h_alpha^3*Nact_p_alpha^2*Nunact_p_alpha + ... % A
                          Nact_h_alpha^2*Nact_h_beta*Nact_p_alpha^2*Nunact_p_beta + ...
                          Nact_h_alpha^2*Nact_h_beta*Nact_p_alpha*Nunact_p_beta*Nact_p_beta + ... % B
                          Nact_h_beta^2*Nact_h_alpha*Nact_p_alpha*Nact_p_beta*Nunact_p_beta + ...
                          Nact_h_alpha*Nact_h_beta^2*Nunact_p_alpha*Nact_p_beta^2 + ... % C
                          Nact_h_beta^3*Nact_p_beta^2*Nunact_p_beta; % D


                % iJK^ABc = iJK^ABc + 
                %           iJK~^AbC~ + iJK~^ABc~  + IJk~^AbC~ +  IJk~^ABc~ +
                %           iJ~K~^aB~C~ + iJ~K~^AB~c~ + Ij~K~^aB~C~ + Ij~K~^AB~c~ +
                %           i~J~K~^A~B~c~
                num_proj4 =  Nunact_h_alpha*Nact_h_alpha*Nact_h_alpha*Nact_p_alpha*Nact_p_alpha*Nunact_p_alpha + ... % A
                             Nunact_h_alpha*Nact_h_alpha*Nact_h_beta*Nact_p_alpha*Nunact_p_beta*Nact_p_beta + ...
                             Nunact_h_alpha*Nact_h_alpha*Nact_h_beta*Nact_p_alpha*Nact_p_alpha*Nunact_p_beta + ...
                             Nact_h_alpha*Nact_h_alpha*Nunact_h_beta*Nact_p_alpha*Nunact_p_alpha*Nact_p_beta + ...
                             Nact_h_alpha*Nact_h_alpha*Nunact_h_beta*Nact_p_alpha*Nact_p_alpha*Nunact_p_beta + ... % B
                             Nunact_h_alpha*Nact_h_beta*Nact_h_beta*Nunact_p_alpha*Nact_p_beta*Nact_p_beta + ...
                             Nunact_h_alpha*Nact_h_beta*Nact_h_beta*Nact_p_alpha*Nact_p_beta*Nunact_p_beta + ...
                             Nact_h_alpha*Nunact_h_beta*Nact_h_beta*Nunact_p_alpha*Nact_p_beta*Nact_p_beta + ...
                             Nact_h_alpha*Nunact_h_beta*Nact_h_beta*Nact_p_alpha*Nact_p_beta*Nunact_p_beta + ... % C
                             Nunact_h_beta*Nact_h_beta*Nact_h_beta*Nact_p_beta*Nact_p_beta*Nunact_p_beta;  % D
            
            sys.size = {[Nvir_alpha, Nocc_alpha], [Nvir_beta, Nocc_beta], ...
                        [Nvir_alpha, Nvir_alpha, Nocc_alpha, Nocc_alpha],...
                        [Nvir_alpha, Nvir_beta, Nocc_alpha, Nocc_beta],...
                        [Nvir_beta, Nvir_beta, Nocc_beta, Nocc_beta],...
                        [Nact_p_alpha, Nact_p_alpha, Nact_p_alpha, Nact_h_alpha, Nact_h_alpha, Nact_h_alpha],...
                        [Nact_p_alpha, Nact_p_alpha, Nunact_p_alpha, Nact_h_alpha, Nact_h_alpha, Nact_h_alpha],...
                        [Nact_p_alpha, Nact_p_alpha, Nact_p_alpha, Nunact_h_alpha, Nact_h_alpha, Nact_h_alpha],...
                        [Nact_p_alpha, Nact_p_alpha, Nunact_p_alpha, Nunact_h_alpha, Nact_h_alpha, Nact_h_alpha],...
                        [Nact_p_alpha, Nact_p_alpha, Nact_p_beta, Nact_h_alpha, Nact_h_alpha, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_alpha, Nunact_p_beta, Nact_h_alpha, Nact_h_alpha, Nact_h_beta],...
                        [Nact_p_alpha, Nunact_p_alpha, Nact_p_beta, Nact_h_alpha, Nact_h_alpha, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_alpha, Nact_p_beta, Nact_h_alpha, Nact_h_alpha, Nunact_h_beta],...
                        [Nact_p_alpha, Nact_p_alpha, Nact_p_beta, Nunact_h_alpha, Nact_h_alpha, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_alpha, Nunact_p_beta, Nunact_h_alpha, Nact_h_alpha, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_alpha, Nunact_p_beta, Nact_h_alpha, Nact_h_alpha, Nunact_h_beta],...
                        [Nact_p_alpha, Nunact_p_alpha, Nact_p_beta, Nunact_h_alpha, Nact_h_alpha, Nact_h_beta],...
                        [Nact_p_alpha, Nunact_p_alpha, Nact_p_beta, Nact_h_alpha, Nact_h_alpha, Nunact_h_beta],...
                        [Nact_p_alpha, Nact_p_beta, Nact_p_beta, Nact_h_alpha, Nact_h_beta, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_beta, Nunact_p_beta, Nact_h_alpha, Nact_h_beta, Nact_h_beta],...
                        [Nunact_p_alpha, Nact_p_beta, Nact_p_beta, Nact_h_alpha, Nact_h_beta, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_beta, Nact_p_beta, Nact_h_alpha, Nunact_h_beta, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_beta, Nact_p_beta, Nunact_h_alpha, Nact_h_beta, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_beta, Nunact_p_beta, Nact_h_alpha, Nunact_h_beta, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_beta, Nunact_p_beta, Nunact_h_alpha, Nact_h_beta, Nact_h_beta],...
                        [Nunact_p_alpha, Nact_p_beta, Nact_p_beta, Nact_h_alpha, Nunact_h_beta, Nact_h_beta],...
                        [Nunact_p_alpha, Nact_p_beta, Nact_p_beta, Nunact_h_alpha, Nact_h_beta, Nact_h_beta],...
                        [Nact_p_beta, Nact_p_beta, Nact_p_beta, Nact_h_beta, Nact_h_beta, Nact_h_beta],...
                        [Nact_p_beta, Nact_p_beta, Nunact_p_beta, Nact_h_beta, Nact_h_beta, Nact_h_beta],...
                        [Nact_p_beta, Nact_p_beta, Nact_p_beta, Nunact_h_beta, Nact_h_beta, Nact_h_beta],...
                        [Nact_p_beta, Nact_p_beta, Nunact_p_beta, Nunact_h_beta, Nact_h_beta, Nact_h_beta]}; 
                    
                                    
           sys.act_triples_dim = Nocc_alpha*Nvir_alpha + ...
                                      Nocc_beta*Nvir_beta + ...
                                      Nocc_alpha^2*Nvir_alpha^2 + ...
                                      Nocc_alpha*Nocc_beta*Nvir_alpha*Nvir_beta + ...
                                      Nocc_beta^2*Nvir_beta^2 + ...
                                      num_proj1 + num_proj2 + num_proj3 + num_proj4;
                     
            if Norb < 40 % allocate active triples 2
                

    
                post3a_PPPHHH = [post2c(end)+1:post2c(end)+Nact_h_alpha^3*Nact_p_alpha^3];
                post3a_PPpHHH = [post3a_PPPHHH(end)+1:post3a_PPPHHH(end)+Nact_h_alpha^3*Nact_p_alpha^2*Nunact_p_alpha];
                post3a_PPPhHH = [post3a_PPpHHH(end)+1:post3a_PPpHHH(end)+Nact_h_alpha^2*Nunact_h_alpha*Nact_p_alpha^3];
                post3a_PPphHH = [post3a_PPPhHH(end)+1:post3a_PPPhHH(end)+Nact_h_alpha^2*Nunact_h_alpha*Nact_p_alpha^2*Nunact_p_alpha];

                post3b_PPPHHH = [post3a_PPphHH(end)+1:post3a_PPphHH(end)+Nact_h_alpha^2*Nact_h_beta*Nact_p_alpha^2*Nact_p_beta];
                post3b_PPpHHH = [post3b_PPPHHH(end)+1:post3b_PPPHHH(end)+Nact_h_alpha^2*Nact_h_beta*Nact_p_alpha^2*Nunact_p_beta];
                post3b_PpPHHH = [post3b_PPpHHH(end)+1:post3b_PPpHHH(end)+Nact_h_alpha^2*Nact_h_beta*Nact_p_alpha*Nunact_p_alpha*Nact_p_beta];
                post3b_PPPHHh = [post3b_PpPHHH(end)+1:post3b_PpPHHH(end)+Nact_h_alpha^2*Nunact_h_beta*Nact_p_alpha^2*Nact_p_beta];
                post3b_PPPhHH = [post3b_PPPHHh(end)+1:post3b_PPPHHh(end)+Nunact_h_alpha*Nact_h_alpha*Nact_h_beta*Nact_p_alpha^2*Nact_p_beta];
                post3b_PPphHH = [post3b_PPPhHH(end)+1:post3b_PPPhHH(end)+Nunact_h_alpha*Nact_h_alpha*Nact_h_beta*Nact_p_alpha^2*Nunact_p_beta];
                post3b_PPpHHh = [post3b_PPphHH(end)+1:post3b_PPphHH(end)+Nact_h_alpha^2*Nunact_h_beta*Nact_p_alpha^2*Nunact_p_beta];
                post3b_PpPhHH = [post3b_PPpHHh(end)+1:post3b_PPpHHh(end)+Nunact_h_alpha*Nact_h_alpha*Nact_h_beta*Nact_p_alpha*Nunact_p_alpha*Nact_p_beta];
                post3b_PpPHHh = [post3b_PpPhHH(end)+1:post3b_PpPhHH(end)+Nact_h_alpha^2*Nunact_h_beta*Nact_p_alpha*Nunact_p_beta*Nact_p_beta];
                
                post3c_PPPHHH = [post3b_PpPHHh(end)+1:post3b_PpPHHh(end)+Nact_h_alpha*Nact_h_beta^2*Nact_p_alpha*Nact_p_beta^2];
                post3c_PPpHHH = [post3c_PPPHHH(end)+1:post3c_PPPHHH(end)+Nact_h_alpha*Nact_h_beta^2*Nact_p_alpha*Nact_p_beta*Nunact_p_beta];
                post3c_pPPHHH = [post3c_PPpHHH(end)+1:post3c_PPpHHH(end)+Nact_h_alpha*Nact_h_beta^2*Nunact_p_alpha*Nact_p_beta^2];
                post3c_PPPHhH = [post3c_pPPHHH(end)+1:post3c_pPPHHH(end)+Nact_h_alpha*Nact_h_beta*Nunact_h_beta*Nact_p_alpha*Nact_p_beta^2];
                post3c_PPPhHH = [post3c_PPPHhH(end)+1:post3c_PPPHhH(end)+Nunact_h_alpha*Nact_h_beta^2*Nact_p_alpha*Nact_p_beta^2];
                post3c_PPpHhH = [post3c_PPPhHH(end)+1:post3c_PPPhHH(end)+Nact_h_alpha*Nact_h_beta*Nunact_h_beta*Nact_p_alpha*Nact_p_beta*Nunact_p_beta];
                post3c_PPphHH = [post3c_PPpHhH(end)+1:post3c_PPpHhH(end)+Nunact_h_alpha*Nact_h_beta^2*Nact_p_alpha*Nact_p_beta*Nunact_p_beta];
                post3c_pPPHhH = [post3c_PPphHH(end)+1:post3c_PPphHH(end)+Nact_h_alpha*Nact_h_beta*Nunact_h_beta*Nunact_p_alpha*Nact_p_beta^2];
                post3c_pPPhHH = [post3c_pPPHhH(end)+1:post3c_pPPHhH(end)+Nunact_h_alpha*Nact_h_beta^2*Nunact_p_alpha*Nact_p_beta^2];
                
                post3d_PPPHHH = [post3c_pPPhHH(end)+1:post3c_pPPhHH(end)+Nact_h_beta^3*Nact_p_beta^3];
                post3d_PPpHHH = [post3d_PPPHHH(end)+1:post3d_PPPHHH(end)+Nact_h_beta^3*Nact_p_beta^2*Nunact_p_beta];
                post3d_PPPhHH = [post3d_PPpHHH(end)+1:post3d_PPpHHH(end)+Nact_h_beta^2*Nunact_h_beta*Nact_p_beta^3];
                post3d_PPphHH = [post3d_PPPhHH(end)+1:post3d_PPPhHH(end)+Nact_h_beta^2*Nunact_h_beta*Nact_p_beta^2*Nunact_p_beta];
           

                sys.posv = {post1a, post1b, post2a, post2b, post2c,...
                            post3a_PPPHHH, post3a_PPpHHH, post3a_PPPhHH, post3a_PPphHH, ...
                            post3b_PPPHHH, post3b_PPpHHH, post3b_PpPHHH, post3b_PPPHHh, post3b_PPPhHH, post3b_PPphHH, post3b_PPpHHh, post3b_PpPhHH, post3b_PpPHHh, ...
                            post3c_PPPHHH, post3c_PPpHHH, post3c_pPPHHH, post3c_PPPHhH, post3c_PPPhHH, post3c_PPpHhH, post3c_PPphHH, post3c_pPPHhH, post3c_pPPhHH, ... 
                            post3d_PPPHHH, post3d_PPpHHH, post3d_PPPhHH, post3d_PPphHH};

            else
                fprintf('Due to problem size, active triples II space not allocated\n')
                sys.posv = {post1a, post1b, post2a, post2b, post2c};
            end
        
        case 3
            
            % IJK^ABC = IJK^ABC + IJK~^ABC~ + IJ~K~^AB~C~ + I~J~K~^A~B~C~
            n1_A = Nact_h_alpha^3*Nact_p_alpha^3;
            n1_B = Nact_h_alpha^2*Nact_h_beta*Nact_p_alpha^2*Nact_p_beta;
            n1_C = Nact_h_alpha*Nact_h_beta^2*Nact_p_alpha*Nact_p_beta^2;
            n1_D = Nact_h_beta^3*Nact_p_beta^3;

            num_proj1 = n1_A + n1_B + n1_C + n1_D;
            
            sys.size = {[Nvir_alpha, Nocc_alpha], [Nvir_beta, Nocc_beta], ...
                        [Nvir_alpha, Nvir_alpha, Nocc_alpha, Nocc_alpha],...
                        [Nvir_alpha, Nvir_beta, Nocc_alpha, Nocc_beta],...
                        [Nvir_beta, Nvir_beta, Nocc_beta, Nocc_beta],...
                        [Nact_p_alpha, Nact_p_alpha, Nact_p_alpha, Nact_h_alpha, Nact_h_alpha, Nact_h_alpha],...
                        [Nact_p_alpha, Nact_p_alpha, Nact_p_beta, Nact_h_alpha, Nact_h_alpha, Nact_h_beta],...
                        [Nact_p_alpha, Nact_p_beta, Nact_p_beta, Nact_h_alpha, Nact_h_beta, Nact_h_beta],...
                        [Nact_p_beta, Nact_p_beta, Nact_p_beta, Nact_h_beta, Nact_h_beta, Nact_h_beta]};
                    
                                        
            sys.act_triples_dim = Nocc_alpha*Nvir_alpha + ...
                                  Nocc_beta*Nvir_beta + ...
                                  Nocc_alpha^2*Nvir_alpha^2 + ...
                                  Nocc_alpha*Nocc_beta*Nvir_alpha*Nvir_beta + ...
                                  Nocc_beta^2*Nvir_beta^2 + ...
                                  num_proj1;

            if Norb < 40 % allocate active triples 3
                

                    post3a = [post2c(end)+1:post2c(end)+n1_A];

                    post3b = [post3a(end)+1:post3a(end)+n1_B];

                    post3c = [post3b(end)+1:post3b(end)+n1_C];

                    post3d = [post3c(end)+1:post3c(end)+n1_D];

                    sys.posv = {post1a, post1b, post2a, post2b, post2c, post3a, post3b, post3c, post3d}; 

            else
                fprintf('Due to problem size, active triples III space not allocated\n')
                sys.posv = {post1a, post1b, post2a, post2b, post2c};
            end
         
        otherwise
            fprintf('Triples active scheme not supported!\n')
            sys.posv = {post1a, post1b, post2a, post2b, post2c};
    end
    

    % diagonal fock masks
    Zocc_alpha = ones(sys.Nocc_alpha) - eye(sys.Nocc_alpha);
    Zvir_alpha = ones(sys.Nvir_alpha) - eye(sys.Nvir_alpha);
    Zocc_beta = ones(sys.Nocc_beta) - eye(sys.Nocc_beta);
    Zvir_beta = ones(sys.Nvir_beta) - eye(sys.Nvir_beta);
    
    % fA
    sys.fa_oo = FA(iocc_alpha, iocc_alpha);  sys.fa_oo_masked = sys.fa_oo.*Zocc_alpha;
    sys.fa_vv = FA(ivir_alpha, ivir_alpha);  sys.fa_vv_masked = sys.fa_vv.*Zvir_alpha;
    sys.fa_ov = FA(iocc_alpha, ivir_alpha);
    sys.fa_vo = FA(ivir_alpha, iocc_alpha);

    % fB
    sys.fb_oo = FB(iocc_beta, iocc_beta);   sys.fb_oo_masked = sys.fb_oo.*Zocc_beta;
    sys.fb_vv = FB(ivir_beta, ivir_beta);   sys.fb_vv_masked = sys.fb_vv.*Zvir_beta;
    sys.fb_ov = FB(iocc_beta, ivir_beta);
    sys.fb_vo = FB(ivir_beta, iocc_beta);

    % vA
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

    % vB
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

    % vC
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
    
    % active space slicing (added as needed)
    
    % diagonal fock masks
    Zact_HH_alpha = ones(sys.Nact_h_alpha) - eye(sys.Nact_h_alpha);
    Zact_PP_alpha = ones(sys.Nact_p_alpha) - eye(sys.Nact_p_alpha);
    Zact_HH_beta = ones(sys.Nact_h_beta) - eye(sys.Nact_h_beta);
    Zact_PP_beta = ones(sys.Nact_p_beta) - eye(sys.Nact_p_beta);
    Zact_hh_alpha = ones(sys.Nunact_h_alpha) - eye(sys.Nunact_h_alpha);
    Zact_pp_alpha = ones(sys.Nunact_p_alpha) - eye(sys.Nunact_p_alpha);
    Zact_hh_beta = ones(sys.Nunact_h_beta) - eye(sys.Nunact_h_beta);
    Zact_pp_beta = ones(sys.Nunact_p_beta) - eye(sys.Nunact_p_beta);
    
    % fA
    sys.fa_HH = FA(iact_h_alpha, iact_h_alpha);  sys.fa_HH_masked = sys.fa_HH.*Zact_HH_alpha;
    sys.fa_PP = FA(iact_p_alpha, iact_p_alpha);  sys.fa_PP_masked = sys.fa_PP.*Zact_PP_alpha;
    sys.fa_HP = FA(iact_h_alpha, iact_p_alpha);
    sys.fa_PH = FA(iact_p_alpha, iact_h_alpha);
    sys.fa_hh = FA(iunact_h_alpha, iunact_h_alpha);  sys.fa_hh_masked = sys.fa_hh.*Zact_hh_alpha;
    sys.fa_pp = FA(iunact_p_alpha, iunact_p_alpha);  sys.fa_pp_masked = sys.fa_pp.*Zact_pp_alpha;
    sys.fa_hp = FA(iunact_h_alpha, iunact_p_alpha);
    sys.fa_ph = FA(iunact_p_alpha, iunact_h_alpha);
    
    % fB
    sys.fb_HH = FB(iact_h_beta, iact_h_beta);   sys.fb_HH_masked = sys.fb_HH.*Zact_HH_beta;
    sys.fb_PP = FB(iact_p_beta, iact_p_beta);   sys.fb_PP_masked = sys.fb_PP.*Zact_PP_beta;
    sys.fb_HP = FB(iact_h_beta, iact_p_beta);
    sys.fb_PH = FB(iact_p_beta, iact_h_beta);
    sys.fb_hh = FB(iunact_h_beta, iunact_h_beta);   sys.fb_hh_masked = sys.fb_hh.*Zact_hh_beta;
    sys.fb_pp = FB(iunact_p_beta, iunact_p_beta);   sys.fb_pp_masked = sys.fb_pp.*Zact_pp_beta;
    sys.fb_hp = FB(iunact_h_beta, iunact_p_beta);
    sys.fb_ph = FB(iunact_p_beta, iunact_h_beta);
    
    % vA
    sys.vA_HHHH = VA(iact_h_alpha, iact_h_alpha, iact_h_alpha, iact_h_alpha);
    sys.vA_HHHP = VA(iact_h_alpha, iact_h_alpha, iact_h_alpha, iact_p_alpha);
    sys.vA_HHPH = VA(iact_h_alpha, iact_h_alpha, iact_p_alpha, iact_h_alpha);
    sys.vA_HPHH = VA(iact_h_alpha, iact_p_alpha, iact_h_alpha, iact_h_alpha);
    sys.vA_PHHH = VA(iact_p_alpha, iact_h_alpha, iact_h_alpha, iact_h_alpha);
    sys.vA_HHPP = VA(iact_h_alpha, iact_h_alpha, iact_p_alpha, iact_p_alpha);
    sys.vA_HPHP = VA(iact_h_alpha, iact_p_alpha, iact_h_alpha, iact_p_alpha);
    sys.vA_PHHP = VA(iact_p_alpha, iact_h_alpha, iact_h_alpha, iact_p_alpha);
    sys.vA_PHPH = VA(iact_p_alpha, iact_h_alpha, iact_p_alpha, iact_h_alpha);
    sys.vA_PPHH = VA(iact_p_alpha, iact_p_alpha, iact_h_alpha, iact_h_alpha);
    sys.vA_HPPH = VA(iact_h_alpha, iact_p_alpha, iact_p_alpha, iact_h_alpha);
    sys.vA_HPPP = VA(iact_h_alpha, iact_p_alpha, iact_p_alpha, iact_p_alpha);
    sys.vA_PHPP = VA(iact_p_alpha, iact_h_alpha, iact_p_alpha, iact_p_alpha);
    sys.vA_PPHP = VA(iact_p_alpha, iact_p_alpha, iact_h_alpha, iact_p_alpha);
    sys.vA_PPPH = VA(iact_p_alpha, iact_p_alpha, iact_p_alpha, iact_h_alpha);
    sys.vA_PPPP = VA(iact_p_alpha, iact_p_alpha, iact_p_alpha, iact_p_alpha);
    
    % vB
    sys.vB_HHHH = VB(iact_h_alpha, iact_h_beta, iact_h_alpha, iact_h_beta);
    sys.vB_HHHP = VB(iact_h_alpha, iact_h_beta, iact_h_alpha, iact_p_beta);
    sys.vB_HHPH = VB(iact_h_alpha, iact_h_beta, iact_p_alpha, iact_h_beta);
    sys.vB_HPHH = VB(iact_h_alpha, iact_p_beta, iact_h_alpha, iact_h_beta);
    sys.vB_PHHH = VB(iact_p_alpha, iact_h_beta, iact_h_alpha, iact_h_beta);
    sys.vB_HHPP = VB(iact_h_alpha, iact_h_beta, iact_p_alpha, iact_p_beta);
    sys.vB_HPHP = VB(iact_h_alpha, iact_p_beta, iact_h_alpha, iact_p_beta);
    sys.vB_PHHP = VB(iact_p_alpha, iact_h_beta, iact_h_alpha, iact_p_beta);
    sys.vB_PHPH = VB(iact_p_alpha, iact_h_beta, iact_p_alpha, iact_h_beta);
    sys.vB_PPHH = VB(iact_p_alpha, iact_p_beta, iact_h_alpha, iact_h_beta);
    sys.vB_HPPH = VB(iact_h_alpha, iact_p_beta, iact_p_alpha, iact_h_beta);
    sys.vB_HPPP = VB(iact_h_alpha, iact_p_beta, iact_p_alpha, iact_p_beta);
    sys.vB_PHPP = VB(iact_p_alpha, iact_h_beta, iact_p_alpha, iact_p_beta);
    sys.vB_PPHP = VB(iact_p_alpha, iact_p_beta, iact_h_alpha, iact_p_beta);
    sys.vB_PPPH = VB(iact_p_alpha, iact_p_beta, iact_p_alpha, iact_h_beta);
    sys.vB_PPPP = VB(iact_p_alpha, iact_p_beta, iact_p_alpha, iact_p_beta);
    
    % vC
    sys.vC_HHHH = VC(iact_h_beta, iact_h_beta, iact_h_beta, iact_h_beta);
    sys.vC_HHHP = VC(iact_h_beta, iact_h_beta, iact_h_beta, iact_p_beta);
    sys.vC_HHPH = VC(iact_h_beta, iact_h_beta, iact_p_beta, iact_h_beta);
    sys.vC_HPHH = VC(iact_h_beta, iact_p_beta, iact_h_beta, iact_h_beta);
    sys.vC_PHHH = VC(iact_p_beta, iact_h_beta, iact_h_beta, iact_h_beta);
    sys.vC_HHPP = VC(iact_h_beta, iact_h_beta, iact_p_beta, iact_p_beta);
    sys.vC_HPHP = VC(iact_h_beta, iact_p_beta, iact_h_beta, iact_p_beta);
    sys.vC_PHHP = VC(iact_p_beta, iact_h_beta, iact_h_beta, iact_p_beta);
    sys.vC_PHPH = VC(iact_p_beta, iact_h_beta, iact_p_beta, iact_h_beta);
    sys.vC_PPHH = VC(iact_p_beta, iact_p_beta, iact_h_beta, iact_h_beta);
    sys.vC_HPPH = VC(iact_h_beta, iact_p_beta, iact_p_beta, iact_h_beta);
    sys.vC_HPPP = VC(iact_h_beta, iact_p_beta, iact_p_beta, iact_p_beta);
    sys.vC_PHPP = VC(iact_p_beta, iact_h_beta, iact_p_beta, iact_p_beta);
    sys.vC_PPHP = VC(iact_p_beta, iact_p_beta, iact_h_beta, iact_p_beta);
    sys.vC_PPPH = VC(iact_p_beta, iact_p_beta, iact_p_beta, iact_h_beta);
    sys.vC_PPPP = VC(iact_p_beta, iact_p_beta, iact_p_beta, iact_p_beta);
    
    % Goldstone e2int slicing
%     sys.e2int_voov = e2int(ivir_alpha,iocc_alpha,iocc_alpha,ivir_alpha);
%     sys.e2int_vovo = e2int(ivir_alpha,iocc_alpha,ivir_alpha,iocc_alpha);

end