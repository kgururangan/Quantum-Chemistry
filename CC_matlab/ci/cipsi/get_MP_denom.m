function [D_MP] = get_MP_denom(occ,unocc,fock)

    D_MP = sum(diag(fock(unocc,unocc))) - sum(diag(fock(occ,coc)));
    
end