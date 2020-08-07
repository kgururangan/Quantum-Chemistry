function [PT1] = evaluate_PT1(Psi, Phi, sys)

    C = Psi.coefs;
    Dets = Psi.dets;
    E = Psi.energy(1);
    
    for L = 1:length(C)
        [det1p, det2p, occ, unocc, sign, excit_rank] = reorder_max_coincidence(Dets(L,:), Phi);
        [denom] = get_MP_denom(occ,unocc,sys.fock);
        PT1 = PT1 + C(L)*slater_eval(sys.e1int,sys.e2int,det1p,det2p,sign,excit_rank)/(E-denom);
    end
    
end

