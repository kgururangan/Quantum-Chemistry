function [MO] = ao_to_mo(AO,C)
    
    if length(size(AO)) == 2
        T1 = einsum(C,AO,'jb,ij->ib');
        MO = einsum(C,T1,'ia,ib->ab');
    else
        T1 = einsum(C,AO,'ld,ijkl->ijkd');
        T2 = einsum(C,T1,'kc,ijkd->ijcd');
        T3 = einsum(C,T2,'jb,ijcd->ibcd');
        MO = einsum(C,T3,'ia,ibcd->abcd');
    end
end
% def ao_to_mo(AO,C):
%     if len(AO.shape) == 2:
%         T1 = np.einsum('jb,ij->ib',C,AO,optimize=True)
%         T2 = np.einsum('ia,ib->ab',C,T1,optimize=True)
%         return T2
%     else:
%         T1 = np.einsum('ld,ijkl->ijkd',C,AO,optimize=True)
%         T2 = np.einsum('kc,ijkd->ijcd',C,T1,optimize=True)
%         T3 = np.einsum('jb,ijcd->ibcd',C,T2,optimize=True)
%         T4 = np.einsum('ia,ibcd->abcd',C,T3,optimize=True)
%         return T4