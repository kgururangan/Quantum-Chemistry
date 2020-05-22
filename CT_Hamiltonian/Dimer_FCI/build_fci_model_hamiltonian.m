function [H] = build_fci_model_hamiltonian(dets, dets_char, J_interact, J_coul, omega_ge)


    fprintf('Building FCI Model Hamiltonian...\n')
    tic

        ndet = size(dets,1);
        H = zeros(ndet);
        for K = 1:ndet
            for L = K:ndet
                if K == L
                    [~, ~, excit_rank] = reorder_max_coincidence(dets(K,:), dets(1,:));
                    H(K,K) = omega_ge*excit_rank;
                    if dets_char(K) == 2
                        H(K,K) = H(K,K) + J_coul;
                    end
                else
                    if dets_char(K) ~= 0 && dets_char(L) ~= 0
                        H(K,L) = J_interact(dets_char(K),dets_char(L));
                        H(L,K) = H(K,L);
                    end
                end
            end
        end
        
    fprintf('FCI Hamiltonian constructed in %4.1f s\n',toc)
    
    
end




