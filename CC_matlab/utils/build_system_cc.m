function [sys] = build_system_v3(e1int,e2int,Nocc)

    ZM = spatial_to_spinorb(e1int);
    
    VM = spatial_to_spinorb(e2int);
    VM = VM - permute(VM,[1,2,4,3]);
    
    occ = 1:Nocc;
    Nunocc = size(ZM,1) - Nocc;
    unocc = Nocc+1:(Nocc+Nunocc);
    
    FM = zeros(size(VM,1));
    for p = 1:size(VM,1)
        for q = 1:size(VM,1)
            FM(p,q) = ZM(p,q);
            for i = 1:Nocc
                FM(p,q) = FM(p,q) + VM(p,i,q,i);
            end
        end
    end
    
    sys.FM = FM;
    sys.VM = VM;
    sys.occ = occ;
    sys.unocc = unocc;
    sys.Nocc = Nocc;
    sys.Nunocc = Nunocc;
    sys.Nelec = Nocc;
    
end