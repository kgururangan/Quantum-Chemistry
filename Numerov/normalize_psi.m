function [ psi_norm ] = normalize_psi( psi )
    psi_norm = psi./(sqrt(sum(conj(psi).*psi)));
end

