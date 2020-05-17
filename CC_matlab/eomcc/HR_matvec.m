function [SIGMA] = HR_matvec(R,HBar)

    [Nocc, Nunocc] = size(HBar{1}{1,2});
    Nov = Nocc*Nunocc;

    SIGMA = zeros(size(R));
    
    r1 = reshape(R(1:Nov),Nunocc,Nocc);
    r2 = reshape(R(Nov+1:end),Nunocc,Nunocc,Nocc,Nocc);
    [sigma_r1, sigma_r2] = build_HBar_R(r1,r2,HBar);
    SIGMA(1:Nov) = sigma_r1(:);
    SIGMA(Nov+1:end) = sigma_r2(:);

end

