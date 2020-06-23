function [HSS] = build_cis_hamiltonian(FM,VM,occ,unocc)


    Nocc = length(occ); Nunocc = length(unocc); Nov = Nocc*Nunocc;
    HF = [ones(1,Nocc),zeros(1,Nunocc)];
    
    [~, s_exc_idx] = get_single_excits(HF,occ,unocc);
    
    % singles-singles block
    HSS = zeros(Nov);
    for p = 1:Nov
        I = s_exc_idx(p,1);
        A = s_exc_idx(p,2);
        for q = 1:Nov
            J = s_exc_idx(q,1);
            B = s_exc_idx(q,2);
            HSS(p,q) = FM(A,B)*(I==J) - FM(I,J)*(A==B) + VM(A,J,I,B);
        end
    end

end

function [s02] = swapPositions(s0,pos1,pos2)
    
    s02 = s0;
    s02(pos1) = s0(pos2);
    s02(pos2) = s0(pos1);
    
end

function [single_exc, single_exc_idx] = get_single_excits(s0, iocc, iunocc)

    Nocc = length(iocc);
    Nunocc = length(iunocc);
    
    single_exc = zeros(Nocc*Nunocc,length(s0));
    single_exc_idx = zeros(Nocc*Nunocc,2);
    
    ct = 1;
    for i = 1:Nocc
        for a = 1:Nunocc
            single_exc(ct,:) = swapPositions(s0,iocc(i),iunocc(a));
            single_exc_idx(ct,:) = [iocc(i), iunocc(a)];
            ct = ct + 1;
        end
    end
    
end
