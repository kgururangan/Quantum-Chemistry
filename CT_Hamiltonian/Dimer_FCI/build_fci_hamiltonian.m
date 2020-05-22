function [H] = build_fci_hamiltonian(ci_vec,z,v,s)

    addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab/utils'));

    fprintf('Building FCI Hamiltonian...\n')
    tic

        if s ~= eye(size(z,1)) % non-orthogonal orbitals
            v = einsum(s,einsum(s,einsum(s,einsum(s,v,'ap,pqrs->aqrs'),'bq,aqrs->abrs'),'cr,abrs->abcs'),'ds,abcs->abcd');
            z = einsum(s,einsum(s,z,'ap,pq->aq'),'bq,aq->ab');
        end

        ndet = size(ci_vec,1);
        H = zeros(ndet);
        for K = 1:ndet
            for L = K:ndet
                [onebody, twobody] = slater_eval(z,v,ci_vec(K,:),ci_vec(L,:));
                H(K,L) = onebody + twobody;
                H(L,K) = H(K,L);
            end
        end
    fprintf('FCI Hamiltonian constructed in %4.1f s\n',toc)
    
    
end

