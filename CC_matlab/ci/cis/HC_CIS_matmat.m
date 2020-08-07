function [SIGMA] = HC_CIS_matmat(C,sys)

    Nocc = sys.Nocc/2; Nunocc = sys.Nunocc/2;

    SIGMA = zeros(size(C));
    
    for i = 1:size(C,2)
        c1 = reshape(C,Nunocc,Nocc);
        sigma_c1 = build_HC_on_singles(c1,sys);
        SIGMA = sigma_c1(:);
    end

end