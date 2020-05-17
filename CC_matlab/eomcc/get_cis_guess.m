function [omega, B] = get_cis_guess(FM,VM,occ,unocc,nroot,nvec_per_root)

    HSS = build_cis_hamiltonian(FM,VM,occ,unocc);
    
    [V,E] = eig(HSS); omega = diag(E); [omega,idx] = sort(omega,'ascend');
    
    V = V(:,idx); omega = omega(1:nroot);
    
    B = V(:,1:nroot*nvec_per_root);
    
end

