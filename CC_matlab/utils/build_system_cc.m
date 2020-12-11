function [sys] = build_system_cc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,NFZ_core,NFZ_vir,Nact_h,Nact_p)

    Nelec = Nocc_a + Nocc_b;

    if ~exist('NFZ_core','var')
        NFZ_core = 0;
    else
        assert(NFZ_core < (Nelec)/2, 'ERROR MSG: Number of frozen core orbitals exceeds number of occupied spinorbitals!')
    end      
    
    if ~exist('NFZ_vir','var')
        NFZ_vir = 0;
    end

    if ~exist('Nact_h','var')
        Nact_h = 10000;
    end
    
    if ~exist('Nact_p','var')
        Nact_p = 10000;
    end

    % Calculate the scf energy
    hf_energy = 0.0;
    for i = 1:Nocc_a
        hf_energy = hf_energy + e1int(i,i);
    end
    for i = 1:Nocc_b
        hf_energy = hf_energy + e1int(i,i);
    end
    for i = 1:Nocc_a
        for j = 1:Nocc_b
            hf_energy = hf_energy + e2int(i,j,i,j);
        end
    end
    for i = 1:Nocc_a
        for j = 1:Nocc_a
            hf_energy = hf_energy + 0.5*(e2int(i,j,i,j)-e2int(i,j,j,i));
        end
    end
    for i = 1:Nocc_b
        for j = 1:Nocc_b
            hf_energy = hf_energy + 0.5*(e2int(i,j,i,j)-e2int(i,j,j,i));
        end
    end

    sys.Escf = hf_energy + Vnuc;
        

    ZM = spatial_to_spinorb(e1int);
    
    VM = spatial_to_spinorb(e2int);
    VM = VM - permute(VM,[1,2,4,3]);
    
    Norb = size(ZM,1);
    sys.Nvir_a = Norb/2 - Nocc_a;
    sys.Nvir_b = Norb/2 - Nocc_b;
    
    Nunocc = Norb - Nelec - 2*NFZ_vir;
    Nocc = Nelec - 2*NFZ_core;
    
    occ = [(2*NFZ_core+1):Nelec];
    unocc = [occ(end)+1:occ(end)+Nunocc];
    
    FM = zeros(size(VM,1));
    for p = 1:size(VM,1)
        for q = 1:size(VM,1)
            FM(p,q) = ZM(p,q);
            for i = 1:Nelec
                FM(p,q) = FM(p,q) + VM(p,i,q,i);
            end
        end
    end

    
    % active space
    M1 = min(Nocc,2*Nact_h);
    M2 = min(Nunocc,2*Nact_p);

    % used for slicing FM and VM
    iact_h = [Nelec-M1+1:Nelec];
    iact_p = [Nelec+1:Nelec+M2];

    sys.e1int = e1int;
    sys.e2int = e2int;
    
    sys.Norb = Norb;
    sys.FM = FM;
    sys.VM = VM;
    sys.ZM = ZM;
    sys.occ_list = occ;
    sys.unocc_list = unocc;
    sys.act_h_list = iact_h;
    sys.act_p_list = iact_p;
    sys.frzc_list = 1:2*NFZ_core;
    sys.frzv_list = [sys.Norb-NFZ_vir+1:sys.Norb];
    sys.Nocc = Nocc;
    sys.Nunocc = Nunocc;
    sys.Nelec = Nelec;
    sys.Nocc_a = Nocc_a;
    sys.Nocc_b = Nocc_b;
    sys.nfzc = 2*NFZ_core;
    sys.nfzv = 2*NFZ_vir;
    sys.Nact_h = M1;
    sys.Nact_p = M2;
    sys.iact_h = [Nocc-M1+1:Nocc]; % used for slicing t, HBar, etc.
    sys.iact_p = [1:M2]; % used for slicing t, HBar, etc.
    sys.Ncore = sys.Nocc-sys.Nact_h;
    sys.Nvir = sys.Nunocc-sys.Nact_p;
    
    % system sizes and dimensions
    sys.Nov = Nunocc*Nocc;
    sys.Noovv = sys.Nov^2;
    sys.Nooovvv = sys.Nov^3;
    sys.singles_dim = sys.Nov;
    sys.doubles_dim = sys.Nov + sys.Noovv;
    sys.triples_dim = sys.Nov + sys.Noovv + sys.Nooovvv;
    sys.triples3_dim = sys.Nov + sys.Noovv + M1^3*M2^3;
    sys.posv1 = 1:sys.singles_dim;
    sys.posv2 = sys.singles_dim+1:sys.doubles_dim;
    sys.posv3 = sys.doubles_dim+1:sys.triples_dim;
    sys.posv3_act = sys.doubles_dim+1:sys.triples3_dim;
    
    % occupied/unocciuped slicing of FM and VM matrices
    sys.foo = FM(occ,occ);
    sys.fvv = FM(unocc,unocc);
    sys.fov = FM(occ,unocc);
    sys.fvo = FM(unocc,occ);
    sys.foo_masked = sys.foo.*(ones(Nocc)-eye(Nocc));
    sys.fvv_masked = sys.fvv.*(ones(Nunocc)-eye(Nunocc));
    
    sys.Voooo = VM(occ,occ,occ,occ);
    sys.Vooov = VM(occ,occ,occ,unocc);
    sys.Voovo = VM(occ,occ,unocc,occ);
    sys.Vovoo = VM(occ,unocc,occ,occ);
    sys.Vvooo = VM(unocc,occ,occ,occ);
    sys.Voovv = VM(occ,occ,unocc,unocc);
    sys.Vovov = VM(occ,unocc,occ,unocc);
    sys.Vvoov = VM(unocc,occ,occ,unocc);
    sys.Vovvo = VM(occ,unocc,unocc,occ);
    sys.Vvvoo = VM(unocc,unocc,occ,occ);
    sys.Vvovo = VM(unocc,occ,unocc,occ);
    sys.Vovvv = VM(occ,unocc,unocc,unocc);
    sys.Vvovv = VM(unocc,occ,unocc,unocc);
    sys.Vvvov = VM(unocc,unocc,occ,unocc);
    sys.Vvvvo = VM(unocc,unocc,unocc,occ);
    sys.Vvvvv = VM(unocc,unocc,unocc,unocc);
    
    % active space slicing
    sys.fHH = FM(iact_h,iact_h);
    sys.fPP = FM(iact_p,iact_p);
    sys.fHP = FM(iact_h,iact_p);
    sys.fPH = FM(iact_p,iact_h);
    sys.fHH_masked = sys.fHH.*(ones(M1)-eye(M1));
    sys.fPP_masked = sys.fPP.*(ones(M2)-eye(M2));
    
    sys.VHHHH = VM(iact_h,iact_h,iact_h,iact_h);
    sys.VHHHP = VM(iact_h,iact_h,iact_h,iact_p);
    sys.VHHPH = VM(iact_h,iact_h,iact_p,iact_h);
    sys.VHPHH = VM(iact_h,iact_p,iact_h,iact_h);
    sys.VPHHH = VM(iact_p,iact_h,iact_h,iact_h);
    sys.VHHPP = VM(iact_h,iact_h,iact_p,iact_p);
    sys.VHPHP = VM(iact_h,iact_p,iact_h,iact_p);
    sys.VPHHP = VM(iact_p,iact_h,iact_h,iact_p);
    sys.VHPPH = VM(iact_h,iact_p,iact_p,iact_h);
    sys.VPPHH = VM(iact_p,iact_p,iact_h,iact_h);
    sys.VPHPH = VM(iact_p,iact_h,iact_p,iact_h);
    sys.VHPPP = VM(iact_h,iact_p,iact_p,iact_p);
    sys.VPHPP = VM(iact_p,iact_h,iact_p,iact_p);
    sys.VPPHP = VM(iact_p,iact_p,iact_h,iact_p);
    sys.VPPPH = VM(iact_p,iact_p,iact_p,iact_h);
    sys.VPPPP = VM(iact_p,iact_p,iact_p,iact_p);
    
    % tons of other active/core or active/core/occ or active/core/unocc 
    % slicings that are added as needed
    sys.VHHvP = VM(iact_h,iact_h,unocc,iact_p);
    sys.VHoPP = VM(iact_h,occ,iact_p,iact_p);
    
end