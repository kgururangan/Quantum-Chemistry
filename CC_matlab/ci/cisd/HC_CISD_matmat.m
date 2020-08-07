function [SIGMA] = HC_CISD_matmat(C,sys)

    Nocc = sys.Nocc; Nunocc = sys.Nunocc;
    Nov = Nocc*Nunocc;

    SIGMA = zeros(size(C));
    
    for i = 1:size(C,2)
        c1 = reshape(C(1:Nov,i),Nunocc,Nocc);
        c2 = reshape(C(Nov+1:end,i),Nunocc,Nunocc,Nocc,Nocc);
        sigma_c1 = build_CISD_HC_on_singles(c1,c2,sys);
        sigma_c2 = build_CISD_HC_on_doubles(c1,c2,sys);
        SIGMA(1:Nov,i) = sigma_c1(:);
        SIGMA(Nov+1:end,i) = sigma_c2(:);
    end

end