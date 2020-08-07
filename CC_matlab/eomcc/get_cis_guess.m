function [omega, B] = get_cis_guess(sys,nroot,nvec_per_root)

    HSS = build_cis_hamiltonian(sys.FM,sys.VM,sys.occ_list-sys.nfzc,sys.unocc_list-sys.nfzc);
    
    [V,E] = eig(HSS); omega = diag(E); [omega,idx] = sort(omega,'ascend');
    
    V = V(:,idx); omega = omega(1:nroot);
    
    B = V(:,1:nroot*nvec_per_root);
    
end

